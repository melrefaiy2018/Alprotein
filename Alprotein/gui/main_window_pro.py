"""
Professional Main Window Layout for Alprotein GUI

This implements the redesigned professional layout with:
- Left: Workflow panel for step-by-step guidance
- Top-Center: 3D viewer with structure info tabs
- Middle-Center: Parameter configuration panel
- Bottom: Multi-tab results display
"""

import os
import sys
from typing import Optional, Dict, Any
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter,
    QTabWidget, QFileDialog, QMessageBox, QGroupBox,
    QLabel, QPushButton, QProgressBar
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QFont

from .widgets.workflow_panel import WorkflowPanel
from .widgets.parameters_panel import ParametersPanel
from .widgets.enhanced_results_panel import EnhancedResultsPanel
from .widgets.protein_viewer import ProteinViewer

# Import Alprotein core components
from ..core.protein_structure import ProteinStructure
from ..core.pigment_system import PigmentSystem
from ..calculators.site_energy_calculator import SiteEnergyCalculator
from ..calculators.hamiltonian_calculator import HamiltonianCalculator


class CalculationWorker(QThread):
    """Worker thread for running calculations in the background"""

    finished = pyqtSignal(str, bool, str)  # calc_type, success, message
    progress = pyqtSignal(int, str)  # value, message

    def __init__(self, calc_type, pigment_system, calculator=None, hamiltonian_calculator=None, params=None):
        super().__init__()
        self.calc_type = calc_type
        self.pigment_system = pigment_system
        self.calculator = calculator
        self.hamiltonian_calculator = hamiltonian_calculator
        self.params = params or {}
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

            if self.calc_type == "site_energies":
                self.progress.emit(10, "Validating system...")

                # Validate system
                validation = self.pigment_system.validate_system()
                if not validation['system_valid']:
                    self.finished.emit(self.calc_type, False,
                                     f"System validation failed: {validation['warnings']}")
                    return

                self.progress.emit(40, "Calculating site energies...")

                # Calculate site energies
                site_energies = self.calculator.calculate_all_site_energies(
                    self.pigment_system, None
                )

                # Build vacuum energies
                vacuum_energies = {}
                for name, pigment in self.pigment_system.items():
                    vacuum_energies[name] = pigment.get_vacuum_energy()

                self.progress.emit(100, "Site energy calculation complete")

                self.results = {
                    'site_energies': site_energies,
                    'vacuum_energies': vacuum_energies,
                    'validation': validation
                }

                message = f"Calculated site energies for {len(site_energies)} pigments"
                self.finished.emit(self.calc_type, True, message)

            elif self.calc_type == "hamiltonian":
                self.progress.emit(20, "Constructing Hamiltonian matrix...")

                # Construct Hamiltonian
                hamiltonian = self.hamiltonian_calculator.construct_hamiltonian()

                self.progress.emit(100, "Hamiltonian constructed")

                self.results = {
                    'hamiltonian': hamiltonian
                }

                message = f"Constructed {hamiltonian.shape[0]}×{hamiltonian.shape[1]} Hamiltonian matrix"
                self.finished.emit(self.calc_type, True, message)

            elif self.calc_type == "diagonalize":
                self.progress.emit(20, "Diagonalizing Hamiltonian...")

                # Diagonalize
                hamiltonian = self.params.get('hamiltonian')
                if hamiltonian is None:
                    raise ValueError("Hamiltonian matrix not provided")

                eigenvalues, eigenvectors = self.hamiltonian_calculator.diagonalize_hamiltonian(hamiltonian)

                self.progress.emit(100, "Diagonalization complete")

                self.results = {
                    'eigenvalues': eigenvalues,
                    'eigenvectors': eigenvectors
                }

                message = f"Computed {len(eigenvalues)} eigenvalues"
                self.finished.emit(self.calc_type, True, message)

            elif self.calc_type == "spectrum":
                self.progress.emit(20, "Calculating absorption and fluorescence spectra...")

                # Get required parameters
                hamiltonian = self.params.get('hamiltonian')
                domains = self.params.get('domains', {})
                temperature = self.params.get('temperature', 300.0)
                disorder_fwhm = self.params.get('disorder_fwhm', 130.0)
                n_ensemble = self.params.get('n_ensemble', 300)

                if hamiltonian is None:
                    raise ValueError("Hamiltonian matrix not provided")

                # Initialize SpectraCalculator
                from Alprotein.calculators.spectra_calculator import SpectraCalculator

                self.progress.emit(30, "Initializing spectrum calculator...")

                spectra_calc = SpectraCalculator(
                    temperature=temperature,
                    disorder_sigma=disorder_fwhm,  # Converted to sigma internally
                    n_ensemble=n_ensemble,
                    dt=0.1,  # fs
                    t_max=2000.0  # fs
                )

                self.progress.emit(50, "Calculating spectra using Renger method...")

                # Calculate both absorption and fluorescence spectra
                wavelengths_abs, absorption, wavelengths_fl, fluorescence = spectra_calc.calculate_spectrum(
                    hamiltonian=hamiltonian,
                    domains=domains,
                    pigment_system=self.pigment_system,
                    include_vibronic=True,
                    use_ensemble=True
                )

                self.progress.emit(100, "Spectra calculation complete")

                self.results = {
                    'wavelengths_abs': wavelengths_abs,
                    'absorption': absorption,
                    'wavelengths_fl': wavelengths_fl,
                    'fluorescence': fluorescence,
                    'method': 'renger',
                    'temperature': temperature,
                    'disorder_fwhm': disorder_fwhm,
                    'n_ensemble': n_ensemble
                }

                message = "Absorption and fluorescence spectra calculated"
                self.finished.emit(self.calc_type, True, message)

            elif self.calc_type == "cdc_analysis":
                self.progress.emit(20, "Starting CDC analysis...")

                # Run detailed CDC analysis
                total_shifts, detailed_contributions = self.calculator.calculate_detailed_site_energy_contributions(
                    self.pigment_system
                )

                self.progress.emit(100, "CDC analysis complete")

                self.results = {
                    'total_shifts': total_shifts,
                    'detailed_contributions': detailed_contributions
                }

                message = f"CDC analysis completed for {len(detailed_contributions)} pigments"
                self.finished.emit(self.calc_type, True, message)

        except Exception as e:
            import traceback
            traceback.print_exc()
            self.finished.emit(self.calc_type, False, str(e))


class MainWindowPro(QWidget):
    """
    Professional main window widget with enhanced layout
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
        self.site_energy_calculator: Optional[SiteEnergyCalculator] = None
        self.hamiltonian_calculator: Optional[HamiltonianCalculator] = None
        self.results: Dict[str, Any] = {}

        # Calculation worker
        self.calc_worker: Optional[CalculationWorker] = None

        # Setup UI
        self.setup_ui()
        self.connect_signals()

        print("🎨 Professional main window initialized")

    def setup_ui(self):
        """Setup the main UI layout"""
        main_layout = QHBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # Create main horizontal splitter
        main_splitter = QSplitter(Qt.Horizontal)
        main_splitter.setStyleSheet("QSplitter::handle { background-color: #e1e1e1; width: 2px; }")

        # Left panel: Workflow
        self.workflow_panel = WorkflowPanel()
        self.workflow_panel.setMinimumWidth(280)
        self.workflow_panel.setMaximumWidth(350)
        main_splitter.addWidget(self.workflow_panel)

        # Right panel: Vertical splitter for top/middle/bottom
        right_splitter = QSplitter(Qt.Vertical)
        right_splitter.setStyleSheet("QSplitter::handle { background-color: #e1e1e1; height: 2px; }")

        # Top: 3D Viewer with tabs
        self.viewer_tab_widget = QTabWidget()
        self.viewer_tab_widget.setStyleSheet("""
            QTabWidget::pane {
                background-color: #ffffff;
                border: 1px solid #e1e1e1;
            }
            QTabBar::tab {
                color: #111111;
                background-color: #f3f3f3;
                padding: 8px 16px;
                border: 1px solid #e1e1e1;
                border-bottom: none;
                border-top-left-radius: 2px;
                border-top-right-radius: 2px;
            }
            QTabBar::tab:selected {
                background-color: #ffffff;
                color: #111111;
                font-weight: bold;
            }
            QTabBar::tab:hover {
                background-color: #ededed;
                color: #111111;
            }
        """)

        self.protein_viewer = ProteinViewer()
        self.viewer_tab_widget.addTab(self.protein_viewer, "3D Structure")

        # Add structure info widget (placeholder for now)
        self.structure_info_widget = self.create_structure_info_widget()
        self.viewer_tab_widget.addTab(self.structure_info_widget, "Structure Info")

        right_splitter.addWidget(self.viewer_tab_widget)

        # Middle: Parameters panel
        self.parameters_panel = ParametersPanel()
        right_splitter.addWidget(self.parameters_panel)

        # Bottom: Results panel
        self.results_panel = EnhancedResultsPanel()
        right_splitter.addWidget(self.results_panel)

        # Set initial sizes: 40% viewer, 30% parameters, 30% results
        right_splitter.setSizes([400, 300, 300])

        main_splitter.addWidget(right_splitter)

        # Set initial horizontal split: 20% workflow, 80% content
        main_splitter.setSizes([280, 1120])

        main_layout.addWidget(main_splitter)

    def create_structure_info_widget(self):
        """Create structure information widget"""
        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.setContentsMargins(20, 20, 20, 20)

        title = QLabel("Structure Information")
        title.setStyleSheet("font-size: 14px; font-weight: bold; color: #111111;")
        layout.addWidget(title)

        self.structure_info_label = QLabel("No structure loaded")
        self.structure_info_label.setStyleSheet("color: #6b6b6b; font-style: italic;")
        self.structure_info_label.setWordWrap(True)
        layout.addWidget(self.structure_info_label)

        layout.addStretch()

        return widget

    def connect_signals(self):
        """Connect signals between components"""
        # Workflow panel signals
        self.workflow_panel.file_selected.connect(self.load_pdb_file_path)
        self.workflow_panel.configure_pigments.connect(self.configure_pigments)
        self.workflow_panel.calculate_site_energies.connect(self.calculate_site_energies)
        self.workflow_panel.construct_hamiltonian.connect(self.construct_hamiltonian)
        self.workflow_panel.diagonalize_hamiltonian.connect(self.diagonalize_hamiltonian)
        self.workflow_panel.calculate_spectrum.connect(self.calculate_spectrum)

        # Parameters panel signals
        self.parameters_panel.parameters_changed.connect(self.on_parameters_changed)

        # Results panel signals
        self.results_panel.pigment_selected.connect(self.protein_viewer.highlight_pigment)

    def load_pdb_file_path(self, file_path: str):
        """Load PDB file from given path"""
        try:
            self.status_message.emit(f"Loading PDB file: {os.path.basename(file_path)}")

            # Load structure
            self.protein_structure = ProteinStructure.from_file_extended(file_path, os.path.basename(file_path))
            self.protein_structure.pdb_path = file_path

            # Create pigment system
            self.pigment_system = PigmentSystem(self.protein_structure)

            if not self.pigment_system.pigments:
                raise ValueError("No pigments found in the structure")

            # Create calculators
            params = self.parameters_panel.get_parameters()
            self.site_energy_calculator = SiteEnergyCalculator(
                dielectric_constant=params.get('dielectric_cdc', 2.0)
            )

            # Create Hamiltonian calculator with all required parameters
            self.hamiltonian_calculator = HamiltonianCalculator(
                self.pigment_system,
                dielectric_cdc=params.get('dielectric_cdc', 2.0),
                dielectric_tresp=params.get('dielectric_tresp', 1.0),
                f_val=params.get('oscillator_strength', 0.72),
                E_0a=14900,  # Default vacuum energy for CLA (cm^-1)
                E_0b=15674   # Default vacuum energy for CHL (cm^-1)
            )

            # Validate system
            validation = self.pigment_system.validate_system()

            # Update UI
            pigment_types = {name: pig.get_resname() for name, pig in self.pigment_system.items()}
            self.workflow_panel.set_pigment_info(len(self.pigment_system), pigment_types)

            # Update structure info
            self.update_structure_info(validation)

            # Load into 3D viewer
            self.protein_viewer.load_structure(self.protein_structure, self.pigment_system)

            # Clear previous results
            self.results = {}
            self.results_panel.clear_results()

            message = f"{len(self.pigment_system)} pigments, {validation['background_atoms']} background atoms"
            self.structure_loaded.emit(True, message)

        except Exception as e:
            error_msg = f"Failed to load PDB file: {str(e)}"
            self.status_message.emit(error_msg)
            self.structure_loaded.emit(False, error_msg)
            QMessageBox.critical(self, "Load Error", error_msg)

    def update_structure_info(self, validation: Dict[str, Any]):
        """Update structure information display"""
        info = f"Structure: {self.protein_structure.name}\n\n"
        info += f"Pigments: {validation.get('total_pigments', 0)}\n"
        info += f"Background atoms: {validation.get('background_atoms', 0)}\n"

        # Get pigment types info
        pigment_types = validation.get('pigment_types', {})
        if pigment_types:
            info += "Pigment types: " + ", ".join([f"{k}({v})" for k, v in pigment_types.items()]) + "\n"

        info += "\n"

        if validation.get('system_valid'):
            info += "✓ System validation: PASSED\n"

            # Show readiness info
            readiness = validation.get('calculation_readiness', {})
            if readiness:
                info += f"  • Site energy ready: {readiness.get('site_energy_ready', 0)}\n"
                info += f"  • Coupling ready: {readiness.get('coupling_ready', 0)}\n"
                info += f"  • Fully ready: {readiness.get('fully_ready', 0)}\n"
        else:
            info += "✗ System validation: FAILED\n"
            warnings = validation.get('warnings', [])
            if warnings:
                info += f"\nWarnings:\n"
                for warning in warnings:
                    info += f"  • {warning}\n"

        self.structure_info_label.setText(info)
        self.structure_info_label.setStyleSheet("color: #111111; background-color: transparent;")

    def configure_pigments(self):
        """Open pigment configuration dialog"""
        QMessageBox.information(self, "Pigment Configuration",
                              "Pigment configuration dialog coming soon!\n\n"
                              "For now, pigments are auto-detected from the structure.")

    def calculate_site_energies(self):
        """Calculate site energies"""
        if not self.pigment_system or not self.site_energy_calculator:
            QMessageBox.warning(self, "No Structure", "Please load a PDB file first.")
            return

        # Start calculation
        self.calc_worker = CalculationWorker(
            "site_energies",
            self.pigment_system,
            calculator=self.site_energy_calculator
        )
        self.calc_worker.finished.connect(self.on_calculation_finished)
        self.calc_worker.progress.connect(self.on_progress_update)

        self.calculation_started.emit("Site Energies")
        self.calc_worker.start()

    def construct_hamiltonian(self):
        """Construct Hamiltonian matrix"""
        if not self.hamiltonian_calculator:
            QMessageBox.warning(self, "No System", "Please load a structure first.")
            return

        # Start calculation
        self.calc_worker = CalculationWorker(
            "hamiltonian",
            self.pigment_system,
            hamiltonian_calculator=self.hamiltonian_calculator
        )
        self.calc_worker.finished.connect(self.on_calculation_finished)
        self.calc_worker.progress.connect(self.on_progress_update)

        self.calculation_started.emit("Hamiltonian")
        self.calc_worker.start()

    def diagonalize_hamiltonian(self):
        """Diagonalize Hamiltonian matrix"""
        if 'hamiltonian' not in self.results:
            QMessageBox.warning(self, "No Hamiltonian",
                              "Please construct the Hamiltonian first.")
            return

        # Start calculation
        self.calc_worker = CalculationWorker(
            "diagonalize",
            self.pigment_system,
            hamiltonian_calculator=self.hamiltonian_calculator,
            params={'hamiltonian': self.results['hamiltonian']['hamiltonian']}
        )
        self.calc_worker.finished.connect(self.on_calculation_finished)
        self.calc_worker.progress.connect(self.on_progress_update)

        self.calculation_started.emit("Diagonalization")
        self.calc_worker.start()

    def calculate_spectrum(self):
        """Calculate absorption and fluorescence spectra"""
        if 'hamiltonian' not in self.results:
            QMessageBox.warning(self, "No Hamiltonian",
                              "Please construct the Hamiltonian first.")
            return

        # Get hamiltonian
        hamiltonian = self.results['hamiltonian']['hamiltonian']

        # Build domains automatically using coupling cutoff
        domain_cutoff = 20.0  # cm⁻¹ - can be made configurable later
        domains = self.hamiltonian_calculator.build_domains(
            cutoff=domain_cutoff,
            hamiltonian=hamiltonian
        )

        # Start calculation
        self.calc_worker = CalculationWorker(
            "spectrum",
            self.pigment_system,
            hamiltonian_calculator=self.hamiltonian_calculator,
            params={
                'hamiltonian': hamiltonian,
                'domains': domains,
                'temperature': 300.0,  # K - can be made configurable later
                'disorder_fwhm': 130.0,  # cm⁻¹ - can be made configurable later
                'n_ensemble': 300  # can be made configurable later
            }
        )
        self.calc_worker.finished.connect(self.on_calculation_finished)
        self.calc_worker.progress.connect(self.on_progress_update)

        self.calculation_started.emit("Spectrum")
        self.calc_worker.start()

    def on_calculation_finished(self, calc_type: str, success: bool, message: str):
        """Handle calculation completion"""
        if success and self.calc_worker and self.calc_worker.results:
            # Store results
            self.results[calc_type] = self.calc_worker.results

            # Update UI based on calculation type
            if calc_type == "site_energies":
                self.results_panel.display_site_energies(self.calc_worker.results)
                self.protein_viewer.update_energy_visualization(
                    self.calc_worker.results['site_energies']
                )
                self.workflow_panel.set_site_energies_complete()

            elif calc_type == "hamiltonian":
                pigment_labels = list(self.pigment_system.pigments.keys())
                self.results_panel.display_hamiltonian(
                    self.calc_worker.results['hamiltonian'],
                    pigment_labels
                )
                self.workflow_panel.set_hamiltonian_constructed()

            elif calc_type == "diagonalize":
                self.results_panel.display_eigenvalues(
                    self.calc_worker.results['eigenvalues'],
                    self.calc_worker.results['eigenvectors']
                )
                self.workflow_panel.set_hamiltonian_diagonalized()

            elif calc_type == "spectrum":
                # Display both absorption and fluorescence spectra
                self.results_panel.display_spectra(
                    self.calc_worker.results['wavelengths_abs'],
                    self.calc_worker.results['absorption'],
                    self.calc_worker.results['wavelengths_fl'],
                    self.calc_worker.results['fluorescence']
                )
                self.workflow_panel.set_spectrum_complete()

            elif calc_type == "cdc_analysis":
                self.results_panel.display_cdc_analysis(self.calc_worker.results)

        # Emit signal for parent window
        self.calculation_finished.emit(calc_type, success, message)

        # Clean up worker
        if self.calc_worker:
            self.calc_worker.deleteLater()
            self.calc_worker = None

    def on_progress_update(self, percentage: int, message: str):
        """Handle progress updates"""
        self.progress_update.emit(percentage, message)

    def on_parameters_changed(self, params: Dict[str, Any]):
        """Handle parameter changes"""
        # Update calculators with new parameters
        if self.site_energy_calculator:
            self.site_energy_calculator.set_parameters(
                dielectric=params.get('dielectric_cdc', 2.0),
                e0a=14900,  # Default vacuum energy for CLA
                e0b=15674   # Default vacuum energy for CHL
            )

        # Update Hamiltonian calculator if it exists
        if self.hamiltonian_calculator:
            # Note: HamiltonianCalculator parameters are set at initialization
            # If parameters change significantly, may need to recreate the calculator
            pass

        self.status_message.emit("Parameters updated")

    def is_calculating(self) -> bool:
        """Check if a calculation is currently running"""
        return self.calc_worker is not None and self.calc_worker.isRunning()

    def stop_calculation(self):
        """Stop the current calculation"""
        if self.calc_worker and self.calc_worker.isRunning():
            self.calc_worker.stop()
            self.calc_worker.wait(3000)
            if self.calc_worker.isRunning():
                self.calc_worker.terminate()

    def cleanup(self):
        """Cleanup resources before closing"""
        self.stop_calculation()

        if hasattr(self, 'protein_viewer'):
            self.protein_viewer.cleanup()

        print("🧹 Professional main window cleanup completed")


# Add numpy import at the top
import numpy as np
