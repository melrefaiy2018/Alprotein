"""
Main Window Layout for Alprotein GUI
"""

import os
import sys
from typing import Optional, Dict, Any
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QSplitter, 
                            QTabWidget, QFileDialog, QMessageBox, QGroupBox,
                            QLabel, QPushButton, QProgressBar, QTextEdit,
                            QTableWidget, QTableWidgetItem, QHeaderView, QSizePolicy)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer, QPropertyAnimation, QRect
from PyQt5.QtGui import QFont

from .widgets.control_sidebar import ControlSidebar
from .widgets.results_panel import ResultsPanel
from .widgets.protein_viewer import ProteinViewer
from .dialogs.settings_dialog import SettingsDialog
from .dialogs.export_dialog import ExportDialog

# Import Alprotein core components
from ..core.protein_structure import ProteinStructure
from ..core.pigment_system import PigmentSystem
from ..calculators.site_energy_calculator import SiteEnergyCalculator


class CalculationWorker(QThread):
    """
    Worker thread for running calculations in the background
    """
    finished = pyqtSignal(str, bool, str)  # calc_type, success, message
    progress = pyqtSignal(int, str)  # value, message
    
    def __init__(self, calc_type, pigment_system, calculator, output_dir=None):
        super().__init__()
        self.calc_type = calc_type
        self.pigment_system = pigment_system
        self.calculator = calculator
        self.output_dir = output_dir
        self.results = None
        self._stop_requested = False
    
    def stop(self):
        """Request the calculation to stop"""
        self._stop_requested = True
    
    def run(self):
        """Run the calculation in the background"""
        try:
            if self._stop_requested:
                return
            
            self.progress.emit(10, "Validating system...")
            
            if self.calc_type == "site_energies":
                # Validate system first
                validation = self.pigment_system.validate_system()
                if not validation['system_valid']:
                    self.finished.emit(self.calc_type, False, 
                                     f"System validation failed: {validation['warnings']}")
                    return
                
                if self._stop_requested:
                    return
                
                self.progress.emit(30, "Calculating site energies...")
                
                # Calculate site energies
                site_energies = self.calculator.calculate_all_site_energies(
                    self.pigment_system, None
                )
                
                # Build vacuum_energies using the actual E0a/E0b used in calculation
                vacuum_energies = {}
                for name, pigment in self.pigment_system.items():
                    pigment_resname = pigment.get_resname().upper()
                    if pigment_resname == 'CLA' and self.calculator.e0a != 0.0:
                        vacuum_energies[name] = self.calculator.e0a
                    elif pigment_resname == 'CHL' and self.calculator.e0b != 0.0:
                        vacuum_energies[name] = self.calculator.e0b
                    else:
                        vacuum_energies[name] = pigment.get_vacuum_energy()
                
                if self._stop_requested:
                    return
                
                self.progress.emit(100, "Site energy calculation complete")
                
                self.results = {
                    'site_energies': site_energies,
                    'vacuum_energies': vacuum_energies,
                    'validation': validation
                }
                
                message = f"Calculated site energies for {len(site_energies)} pigments"
                self.finished.emit(self.calc_type, True, message)
                
            elif self.calc_type == "cdc_analysis":
                self.progress.emit(20, "Starting CDC analysis...")
                
                try:
                    # Run detailed CDC analysis
                    total_shifts, detailed_contributions = self.calculator.calculate_detailed_site_energy_contributions(
                        self.pigment_system
                    )
                    
                    if self._stop_requested:
                        return
                    
                    self.progress.emit(80, "Generating CDC reports...")
                    
                    # If output directory is provided, save detailed analysis
                    if self.output_dir:
                        pdb_path = getattr(self.pigment_system.structure, 'pdb_path', None)
                        if pdb_path:
                            try:
                                detailed_results = self.calculator.perform_cdc_analysis(
                                    self.pigment_system, pdb_path, self.output_dir
                                )
                                self.results = {
                                    'total_shifts': total_shifts,
                                    'detailed_contributions': detailed_contributions,
                                    'detailed_results': detailed_results,
                                    'output_dir': self.output_dir
                                }
                            except Exception as cdc_error:
                                # If detailed analysis fails, still return the basic results
                                print(f"Warning: Could not generate detailed CDC reports: {cdc_error}")
                                self.results = {
                                    'total_shifts': total_shifts,
                                    'detailed_contributions': detailed_contributions
                                }
                        else:
                            self.results = {
                                'total_shifts': total_shifts,
                                'detailed_contributions': detailed_contributions
                            }
                    else:
                        self.results = {
                            'total_shifts': total_shifts,
                            'detailed_contributions': detailed_contributions
                        }
                    
                    self.progress.emit(100, "CDC analysis complete")
                    
                    message = f"CDC analysis completed for {len(detailed_contributions)} pigments"
                    self.finished.emit(self.calc_type, True, message)
                    
                except Exception as cdc_error:
                    error_msg = f"CDC analysis failed: {str(cdc_error)}"
                    print(f"CDC Error Details: {cdc_error}")
                    import traceback
                    traceback.print_exc()
                    self.finished.emit(self.calc_type, False, error_msg)
                    return
                
        except Exception as e:
            self.finished.emit(self.calc_type, False, str(e))


class MainWindow(QWidget):
    """
    Main window widget containing all GUI components
    """
    
    # Signals for communication with parent
    structure_loaded = pyqtSignal(bool, str)  # success, message
    calculation_started = pyqtSignal(str)  # calc_type
    calculation_finished = pyqtSignal(str, bool, str)  # calc_type, success, message
    status_message = pyqtSignal(str)  # message
    progress_update = pyqtSignal(int, str)  # value, message
    
    def __init__(self):
        super().__init__()
        
        # Initialize data
        self.protein_structure: Optional[ProteinStructure] = None
        self.pigment_system: Optional[PigmentSystem] = None
        self.calculator: Optional[SiteEnergyCalculator] = None
        self.results: Dict[str, Any] = {}
        
        # Calculation worker
        self.calc_worker: Optional[CalculationWorker] = None
        
        # Setup UI
        self.setWindowTitle("Alprotein - Protein Site Energy Calculator")
        self.setup_ui()
        self.connect_signals()
        
        print("📋 Main window initialized")
    
    def setup_ui(self):
        """Setup the main UI layout based on the MVP redesign plan."""
        # Create main splitter
        main_splitter = QSplitter(Qt.Horizontal)
        main_splitter.setStyleSheet("QSplitter::handle { background-color: #e1e1e1; }") # Optional styling

        # Left panel: ControlSidebar
        self.control_sidebar = self.create_control_panel()
        main_splitter.addWidget(self.control_sidebar)

        # Right panel: QTabWidget for visualization and results
        right_panel = self.create_visualization_panel()
        main_splitter.addWidget(right_panel)

        # Set initial sizes for the splitter (sidebar, main_panel)
        main_splitter.setSizes([300, 1200])

        # Set the splitter as the main layout
        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.addWidget(main_splitter)

        # Placeholder for floating toolbar
        self.init_floating_toolbar()

    def init_floating_toolbar(self):
        """Initializes a placeholder for a floating toolbar over the 3D viewer."""
        # This method is a placeholder for future implementation.
        # For now, it can be empty or log a message.
        print("🔧 Floating toolbar initialized (placeholder)")
    
    def create_control_panel(self):
        """Create the left control panel"""
        self.control_sidebar = ControlSidebar()
        return self.control_sidebar
    
    def create_visualization_panel(self):
        """Create the right visualization panel"""
        # Create tab widget for different views
        tab_widget = QTabWidget()
        
        # Set tab widget styling to ensure black text
        tab_widget.setStyleSheet("""
            QTabWidget::pane {
                background-color: white;
                border: 1px solid #ddd;
            }
            QTabBar::tab {
                color: black !important;
                background-color: #f0f0f0;
                padding: 8px 16px;
                margin: 2px;
                border: 1px solid #ddd;
                border-bottom: none;
            }
            QTabBar::tab:selected {
                background-color: white;
                color: black !important;
                border-bottom: 2px solid #3399ff;
            }
            QTabBar::tab:hover {
                background-color: #e0e0e0;
                color: black !important;
            }
        """)
        
        # 3D Visualization tab
        self.protein_viewer = ProteinViewer()
        tab_widget.addTab(self.protein_viewer, "3D Structure")
        
        # Results tab
        self.results_panel = ResultsPanel()
        tab_widget.addTab(self.results_panel, "Results & Analysis")
        
        return tab_widget
    
    def connect_signals(self):
        """Connect signals between components"""
        # Connect signals from the consolidated control sidebar
        self.control_sidebar.file_selected.connect(self.load_pdb_file_path)
        self.control_sidebar.calculate_site_energies.connect(self.calculate_site_energies)
        self.control_sidebar.run_cdc_analysis.connect(self.run_cdc_analysis)
        self.control_sidebar.settings_changed.connect(self.update_calculator_settings)
        self.control_sidebar.view_changed.connect(self.protein_viewer.viewer_3d.update_view_settings)
        
        # Results panel signals
        self.results_panel.pigment_selected.connect(self.protein_viewer.highlight_pigment)
        self.results_panel.export_requested.connect(self.export_results)
    
    def load_pdb_file(self):
        """Open file dialog to load PDB file"""
        # Use non-native dialog on macOS to avoid warnings
        options = QFileDialog.DontUseNativeDialog if sys.platform == 'darwin' else QFileDialog.Options()
        
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Load PDB File",
            "",
            "PDB Files (*.pdb);;All Files (*)",
            options=options
        )
        
        if file_path:
            self.load_pdb_file_path(file_path)
    
    def load_pdb_file_path(self, file_path: str):
        """Load PDB file from given path"""
        try:
            self.status_message.emit(f"Loading PDB file: {os.path.basename(file_path)}")
            
            # Load structure using existing method
            self.protein_structure = ProteinStructure.from_file_extended(file_path, os.path.basename(file_path))
            
            # Store the file path for later use
            self.protein_structure.pdb_path = file_path
            
            # Create pigment system
            self.pigment_system = PigmentSystem(self.protein_structure)
            
            # Check if we have pigments
            if not self.pigment_system.pigments:
                raise ValueError("No pigments found in the structure")
            
            # Create calculator with current settings from panel
            current_settings = self.control_sidebar.get_calculation_settings()
            self.calculator = SiteEnergyCalculator(
                dielectric_constant=current_settings.get('dielectric_constant', 2.5),
                e0a=current_settings.get('e0a', 0.0),
                e0b=current_settings.get('e0b', 0.0)
            )
            
            # Validate the system
            validation = self.pigment_system.validate_system()
            
            # Update UI components
            self.control_sidebar.set_loaded_file(file_path, len(self.pigment_system))
            self.control_sidebar.set_system_info(validation)
            self.control_sidebar.update_pigment_list(self.pigment_system.pigments)
            self.protein_viewer.load_structure(self.protein_structure, self.pigment_system)
            
            # Clear previous results
            self.results = {}
            self.results_panel.clear_results()
            
            # Reset calculation panel state
            self.control_sidebar.set_site_energies_calculated(False)
            
            message = f"{len(self.pigment_system)} pigments, {validation['background_atoms']} background atoms"
            self.structure_loaded.emit(True, message)
            
        except Exception as e:
            error_msg = f"Failed to load PDB file: {str(e)}"
            self.status_message.emit(error_msg)
            self.structure_loaded.emit(False, error_msg)
    
    def calculate_site_energies(self):
        """Run site energy calculation"""
        if not self.pigment_system or not self.calculator:
            QMessageBox.warning(self, "No Structure", "Please load a PDB file first.")
            return
        
        # Start calculation in background thread
        self.calc_worker = CalculationWorker(
            "site_energies", self.pigment_system, self.calculator
        )
        self.calc_worker.finished.connect(self.on_calculation_finished)
        self.calc_worker.progress.connect(self.on_progress_update)
        
        self.control_sidebar.set_calculation_running(True, "site_energies")
        self.calculation_started.emit("Site Energies")
        self.calc_worker.start()
    
    def run_cdc_analysis(self):
        """Run CDC analysis"""
        if not self.pigment_system or not self.calculator:
            QMessageBox.warning(self, "No Structure", "Please load a PDB file first.")
            return
        
        # Check if site energies have been calculated
        if 'site_energies' not in self.results:
            reply = QMessageBox.question(
                self, 
                "No Site Energies", 
                "Site energies must be calculated first. Would you like to calculate them now?",
                QMessageBox.Yes | QMessageBox.No
            )
            if reply == QMessageBox.Yes:
                self.calculate_site_energies()
            return
        
        # Ask if user wants detailed file output (optional)
        reply = QMessageBox.question(
            self,
            "CDC Analysis Options",
            "Do you want to generate detailed file reports?\n\n"
            "• Yes: Select output directory for detailed analysis files\n"
            "• No: Run analysis for visualization only",
            QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel
        )
        
        if reply == QMessageBox.Cancel:
            return
        
        output_dir = None
        if reply == QMessageBox.Yes:
            # Ask user for output directory only if they want file output
            output_dir = QFileDialog.getExistingDirectory(
                self, "Select Output Directory for CDC Analysis Reports",
                options=QFileDialog.ShowDirsOnly | QFileDialog.DontUseNativeDialog
            )
            
            if not output_dir:
                return  # User cancelled
        
        # Start calculation in background thread
        self.calc_worker = CalculationWorker(
            "cdc_analysis", self.pigment_system, self.calculator, output_dir
        )
        self.calc_worker.finished.connect(self.on_calculation_finished)
        self.calc_worker.progress.connect(self.on_progress_update)
        
        self.control_sidebar.set_calculation_running(True, "cdc_analysis")
        self.calculation_started.emit("CDC Analysis")
        self.calc_worker.start()
    
    def on_calculation_finished(self, calc_type: str, success: bool, message: str):
        """Handle calculation completion"""
        if success and self.calc_worker and self.calc_worker.results:
            # Store results
            self.results[calc_type] = self.calc_worker.results
            
            # Update results panel
            if calc_type == "site_energies":
                # Include calculator parameters in results for proper display
                self.calc_worker.results['calculator_params'] = {
                    'dielectric_constant': self.calculator.dielectric_constant,
                    'e0a': self.calculator.e0a,
                    'e0b': self.calculator.e0b
                }
                
                self.results_panel.display_site_energies(self.calc_worker.results)
                # Update 3D visualization with energies
                self.protein_viewer.update_energy_visualization(
                    self.calc_worker.results['site_energies']
                )
                # Update action buttons to enable CDC analysis
                self.control_sidebar.set_site_energies_calculated(True)
                
            elif calc_type == "cdc_analysis":
                self.results_panel.display_cdc_analysis(self.calc_worker.results)
                # Update 3D visualization with contributions
                self.protein_viewer.update_contribution_visualization(
                    self.calc_worker.results['detailed_contributions']
                )
        
        # Emit signal for parent window
        self.calculation_finished.emit(calc_type, success, message)
        
        # Reset calculation running state
        self.control_sidebar.set_calculation_running(False)
        
        # Clean up worker
        if self.calc_worker:
            self.calc_worker.deleteLater()
            self.calc_worker = None
    
    def on_progress_update(self, percentage: int, message: str):
        """Handle progress updates from calculation worker"""
        self.control_sidebar.update_progress(percentage, message)
        self.progress_update.emit(percentage, message)
    
    def update_calculator_settings(self, settings: Dict[str, Any]):
        """Update calculator settings"""
        if self.calculator:
            self.calculator.set_parameters(
                dielectric=settings.get('dielectric_constant', 2.5),
                e0a=settings.get('e0a', 0.0),
                e0b=settings.get('e0b', 0.0)
            )
            self.status_message.emit(f"Updated calculator settings: ε={settings.get('dielectric_constant', 2.5):.3f}, E₀a={settings.get('e0a', 0.0):.1f}, E₀b={settings.get('e0b', 0.0):.1f}")

    
    
    def export_results(self):
        """Export calculation results"""
        if not self.results:
            QMessageBox.warning(self, "No Results", "No calculation results to export.")
            return
        
        # Open export dialog
        dialog = ExportDialog(self.results, self)
        if dialog.exec_() == ExportDialog.Accepted:
            export_settings = dialog.get_export_settings()
            
            # Perform export based on settings
            try:
                self.perform_export(export_settings)
                print(f"[INFO] Results exported to: {export_settings['output_dir']}")
            except Exception as e:
                QMessageBox.critical(self, "Export Error", f"Failed to export results:\n{str(e)}")
    
    def perform_export(self, export_settings: Dict[str, Any]):
        """Perform the actual export operation"""
        import json
        import csv
        import numpy as np
        from pathlib import Path
        
        def convert_to_serializable(obj):
            """Convert numpy types and other non-serializable objects to Python types"""
            if isinstance(obj, np.bool_):
                return bool(obj)
            elif isinstance(obj, (np.integer, np.int64, np.int32)):
                return int(obj)
            elif isinstance(obj, (np.floating, np.float64, np.float32)):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, dict):
                return {k: convert_to_serializable(v) for k, v in obj.items()}
            elif isinstance(obj, (list, tuple)):
                return [convert_to_serializable(item) for item in obj]
            else:
                return obj
        
        output_dir = Path(export_settings['output_dir'])
        output_dir.mkdir(exist_ok=True)
        
        # Export site energies if available
        if 'site_energies' in self.results and export_settings.get('export_site_energies', True):
            site_energies = self.results['site_energies']['site_energies']
            
            # CSV export
            if export_settings.get('format_csv', True):
                csv_path = output_dir / 'site_energies.csv'
                with open(csv_path, 'w', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow(['Pigment', 'Site Energy (cm-1)', 'Vacuum Energy (cm-1)', 'Energy Shift (cm-1)'])
                    
                    vacuum_energies = self.results['site_energies']['vacuum_energies']
                    for pigment, site_energy in site_energies.items():
                        vacuum_energy = vacuum_energies.get(pigment, 0.0)
                        energy_shift = site_energy - vacuum_energy
                        writer.writerow([pigment, f"{site_energy:.1f}", f"{vacuum_energy:.1f}", f"{energy_shift:.1f}"])
            
            # JSON export
            if export_settings.get('format_json', True):
                json_path = output_dir / 'site_energies.json'
                with open(json_path, 'w') as f:
                    # Convert to serializable format before saving
                    serializable_data = convert_to_serializable(self.results['site_energies'])
                    json.dump(serializable_data, f, indent=2)
        
        # Export CDC analysis if available
        if 'cdc_analysis' in self.results and export_settings.get('export_cdc', True):
            cdc_data = self.results['cdc_analysis']
            
            # JSON export
            if export_settings.get('format_json', True):
                json_path = output_dir / 'cdc_analysis.json'
                with open(json_path, 'w') as f:
                    # Convert to serializable format before saving
                    serializable_data = convert_to_serializable(cdc_data['detailed_contributions'])
                    json.dump(serializable_data, f, indent=2)
        
        self.status_message.emit(f"Results exported to {output_dir}")
    
    
    
    def is_calculating(self) -> bool:
        """Check if a calculation is currently running"""
        return self.calc_worker is not None and self.calc_worker.isRunning()
    
    def stop_calculation(self):
        """Stop the current calculation"""
        if self.calc_worker and self.calc_worker.isRunning():
            self.calc_worker.stop()
            self.calc_worker.wait(3000)  # Wait up to 3 seconds
            if self.calc_worker.isRunning():
                self.calc_worker.terminate()
    
    def cleanup(self):
        """Cleanup resources before closing"""
        self.stop_calculation()
        
        # Clean up 3D viewer
        if hasattr(self, 'protein_viewer'):
            self.protein_viewer.cleanup()
        
        print("🧹 Main window cleanup completed")
