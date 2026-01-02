"""
Scientific Workbench Window - Complete Integration

Integrates all workbench components into a professional scientific interface:
- ProteinViewer (3D visualization)
- AnalysisDashboard (mini-plots)
- DataTableWidget (data table)
- ToolsPanel (settings and controls)
"""

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter,
    QMessageBox, QFileDialog, QTabWidget, QLabel, QPushButton,
    QDoubleSpinBox, QGroupBox, QFrame, QGridLayout, QSizePolicy,
    QScrollArea, QCheckBox
)
from PyQt5.QtCore import Qt, pyqtSignal, QThread
from PyQt5.QtGui import QFont
import numpy as np
import logging
from pathlib import Path
from typing import Dict, Any, Optional, List

logger = logging.getLogger(__name__)

from Alprotein.gui.widgets.tools_panel import ToolsPanel
from Alprotein.gui.widgets.analysis_dashboard import AnalysisDashboard, MiniPlotWidget as MiniPlot
from Alprotein.gui.widgets.data_table_widget import DataTableWidget
from Alprotein.gui.widgets.protein_viewer import ProteinViewer
from Alprotein.gui.widgets.hamiltonian_widget import HamiltonianWidget
from Alprotein.gui.widgets.spectrum_widget import SpectrumPlotWidget
# from Alprotein.gui.widgets.summary_cards_widget import SummaryCardsWidget
from Alprotein.core.pigment_system import PigmentSystem
from Alprotein.core.protein_structure import ProteinStructure
from Alprotein.calculators.site_energy_calculator import SiteEnergyCalculator
from Alprotein.calculators.hamiltonian_calculator import HamiltonianCalculator
from Alprotein.calculators.exciton_calculator import ExcitonCalculator


# PRD Color Palette
COLORS = {
    'accent_primary': '#111111',      # Primary ink
    'accent_secondary': '#3a3a3a',    # Secondary ink
    'bg_app': '#f7f7f7',              # Light gray
    'bg_panel': '#ffffff',            # Pure white
    'bg_sidebar': '#f3f3f3',          # Soft gray
    'text_primary': '#111111',        # Near-black
    'text_secondary': '#6b6b6b',      # Medium gray
    'text_disabled': '#b0b0b0',       # Light gray
    'success': '#1f7a44',             # Muted green
    'warning': '#8a6d1f',             # Muted amber
    'error': '#8b2b2b',               # Muted red
    'border': '#e1e1e1',              # Border color
}


class CalculationWorker(QThread):
    """Worker thread for background calculations"""

    finished = pyqtSignal(str, object)  # calc_type, result
    error = pyqtSignal(str, str)  # calc_type, error_message
    progress = pyqtSignal(str, int)  # message, percentage

    def __init__(self, calc_type, calculator, *args, **kwargs):
        super().__init__()
        self.calc_type = calc_type
        self.calculator = calculator
        self.args = args
        self.kwargs = kwargs

    def run(self):
        """Run the calculation"""
        try:
            self.progress.emit(f"Running {self.calc_type}...", 10)

            if self.calc_type == "site_energies":
                # Calculate site energies for all pigments
                pigment_system = self.kwargs.get('pigment_system')
                vacuum_energies = {'CLA': 14900.0, 'CHL': 15674.0}
                result = self.calculator.calculate_all_site_energies(
                    pigment_system,
                    vacuum_energies
                )
                self.progress.emit("Site energies calculated", 100)

            elif self.calc_type == "hamiltonian":
                try:
                    logger.info("Starting Hamiltonian construction...")
                    # Pass parameters for coupling calculation
                    params = {'coupling_calc': self.kwargs.get('coupling_calc', 'tresp')}
                    result = self.calculator.construct_hamiltonian(params=params)
                    logger.info(f"Hamiltonian constructed successfully. Type: {type(result)}, Shape: {result.shape}")
                    self.progress.emit("Hamiltonian constructed", 100)
                except Exception as e:
                    import traceback
                    error_details = traceback.format_exc()
                    logger.error(f"Hamiltonian construction failed: {type(e).__name__}: {str(e)}")
                    logger.error(f"Full traceback:\n{error_details}")
                    raise RuntimeError(
                        f"Failed to construct Hamiltonian: {str(e)}\n\n"
                        f"Details:\n{error_details}"
                    ) from e

            elif self.calc_type == "spectrum_renger":
                # Advanced Renger lineshape calculation
                self.progress.emit("Calculating absorption and fluorescence...", 10)

                hamiltonian = self.kwargs['hamiltonian']
                domains = self.kwargs['domains']
                pigment_system = self.kwargs['pigment_system']
                include_vibronic = self.kwargs.get('include_vibronic', True)
                use_ensemble = self.kwargs.get('use_ensemble', True)
                progress_callback = self.kwargs.get('progress_callback')
                if progress_callback is None:
                    def progress_callback(message, percentage):
                        self.progress.emit(message, percentage)

                # Calculate both absorption and fluorescence spectra using Renger theory
                wavelengths_abs, absorption, wavelengths_fl, fluorescence = self.calculator.calculate_spectrum(
                    hamiltonian=hamiltonian,
                    domains=domains,
                    pigment_system=pigment_system,
                    include_vibronic=include_vibronic,
                    use_ensemble=use_ensemble,
                    progress_callback=progress_callback
                )

                self.progress.emit("Spectra calculated (Renger method)", 100)

                result = {
                    'wavelengths_abs': wavelengths_abs,
                    'absorption': absorption,
                    'wavelengths_fl': wavelengths_fl,
                    'fluorescence': fluorescence,
                    'method': 'renger',
                    # For backward compatibility, also include old format (using absorption data)
                    'wavelengths': wavelengths_abs,
                    'spectrum': absorption
                }

            elif self.calc_type == "spectrum":
                # Simple spectrum calculation (fallback/legacy)
                eigenvalues = self.kwargs.get('eigenvalues')
                temperature = self.kwargs.get('temperature', 300)
                broadening = self.kwargs.get('broadening', 20.0)  # nm

                # Convert eigenvalues to wavelengths (cm^-1 to nm)
                wavelengths_nm = 1e7 / eigenvalues  # stick spectrum positions

                # Create broadened spectrum
                wl_range = np.linspace(wavelengths_nm.min() - 100, wavelengths_nm.max() + 100, 1000)
                spectrum = np.zeros_like(wl_range)

                # Equal intensities for now (could be weighted by oscillator strengths)
                for wl in wavelengths_nm:
                    # Gaussian broadening
                    sigma = broadening / (2 * np.sqrt(2 * np.log(2)))
                    spectrum += np.exp(-((wl_range - wl) ** 2) / (2 * sigma ** 2))

                # Normalize
                spectrum /= spectrum.max()

                result = {
                    'wavelengths': wl_range,
                    'spectrum': spectrum,
                    'stick_wavelengths': wavelengths_nm,
                    'stick_intensities': np.ones_like(wavelengths_nm),
                    'method': 'simple'
                }
                self.progress.emit("Spectrum calculated", 100)

            elif self.calc_type == "exciton_distribution":
                # Calculate exciton distributions
                self.progress.emit("Calculating exciton distributions...", 10)

                hamiltonian = self.kwargs['hamiltonian']
                domains = self.kwargs['domains']
                pigment_system = self.kwargs['pigment_system']

                # Calculate distributions
                distributions, labels = self.calculator.calculate_distributions(
                    hamiltonian=hamiltonian,
                    domains=domains,
                    pigment_system=pigment_system
                )

                result = {
                    'distributions': distributions,
                    'labels': labels
                }

                self.progress.emit("Exciton distributions calculated", 100)

            else:
                raise ValueError(f"Unknown calculation type: {self.calc_type}")

            self.finished.emit(self.calc_type, result)

        except Exception as e:
            import traceback
            error_msg = f"{type(e).__name__}: {str(e)}\n\nTraceback:\n{traceback.format_exc()}"
            logger.error(f"Worker thread error in {self.calc_type}: {error_msg}")
            self.error.emit(self.calc_type, error_msg)


class ScientificWorkbenchWindow(QWidget):
    """
    Complete Scientific Workbench integrating all components
    """

    # Signals
    status_message = pyqtSignal(str)
    progress_update = pyqtSignal(str, int)

    def __init__(self, parent=None):
        super().__init__(parent)

        # Core data
        self.pigment_system: Optional[PigmentSystem] = None
        self.site_energy_calculator: Optional[SiteEnergyCalculator] = None
        self.hamiltonian_calculator: Optional[HamiltonianCalculator] = None
        self.exciton_calculator: Optional[ExcitonCalculator] = None

        # Results storage
        self.site_energies: Dict[str, float] = {}
        self.calculated_site_energies: Dict[str, float] = {}
        self.site_energy_overrides: Dict[str, float] = {}
        self.vacuum_energies: Dict[str, float] = {}
        self.hamiltonian: Optional[np.ndarray] = None
        self.eigenvalues: Optional[np.ndarray] = None
        self.eigenvectors: Optional[np.ndarray] = None
        self.spectrum_data: Optional[Dict[str, np.ndarray]] = None
        self.exciton_distributions: Optional[Dict] = None
        self.exciton_labels: Optional[List[str]] = None
        self.pigment_labels: Optional[List[str]] = None

        # Axes preferences for exciton plots
        self.wavelength_min_nm: float = 600.0
        self.wavelength_max_nm: float = 720.0

        # Worker threads
        self.current_worker: Optional[CalculationWorker] = None

        # Workflow control
        self.run_all_flag = False  # Track if running complete workflow

        # Current file
        self.current_file: Optional[Path] = None

        self.setup_ui()
        self.connect_signals()

    def setup_ui(self):
        """Setup the UI layout"""
        main_layout = QHBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # Main horizontal splitter (workspace | tools)
        main_splitter = QSplitter(Qt.Horizontal)

        # LEFT: Workspace with tabs (3D viewer, dashboard, table)
        self.workspace_tabs = QTabWidget()
        # Stylesheet removed to use global styles

        # Tab 0: 3D Structure Viewer (KEEP)
        self.protein_viewer = ProteinViewer()
        self.workspace_tabs.addTab(self.protein_viewer, "🔬 3D Structure")

        # Tab 1: Hamiltonian (NEW)
        self.hamiltonian_widget = HamiltonianWidget()
        self.workspace_tabs.addTab(self.hamiltonian_widget, "📐 Hamiltonian")

        # Tab 2: Spectra (NEW)
        self.spectrum_widget = SpectrumPlotWidget()
        self.workspace_tabs.addTab(self.spectrum_widget, "📊 Spectra")

        # Tab 3: Data Analysis (ENHANCED)
        self.data_analysis_tab = self.create_data_analysis_tab()
        self.workspace_tabs.addTab(self.data_analysis_tab, "📋 Data Analysis")

        # LEFT: Tools Panel (sidebar) - MOVED TO LEFT per PRD
        self.tools_panel = ToolsPanel()
        main_splitter.addWidget(self.tools_panel)

        # RIGHT: Workspace tabs - MOVED TO RIGHT per PRD
        main_splitter.addWidget(self.workspace_tabs)

        # Set proportions: 20% tools (left), 80% workspace (right)
        main_splitter.setSizes([300, 1200])
        main_splitter.setStretchFactor(0, 1)  # Tools panel (left) - lower stretch
        main_splitter.setStretchFactor(1, 4)  # Workspace (right) - higher stretch

        main_layout.addWidget(main_splitter)

    def create_data_analysis_tab(self):
        """Create enhanced Data Analysis tab with professional dashboard layout"""
        tab_widget = QWidget()
        tab_widget.setStyleSheet("background-color: #f7f7f7;")  # Light neutral background
        
        # Main grid layout for the tab content
        main_layout = QVBoxLayout(tab_widget)
        main_layout.setContentsMargins(20, 20, 20, 20)
        main_layout.setSpacing(20)

        # 1. Header Row
        header_container = QWidget()
        header_layout = QHBoxLayout(header_container)
        header_layout.setContentsMargins(0, 0, 0, 0)
        
        # Title
        title_label = QLabel("Exciton Distribution Analysis")
        title_label.setStyleSheet("font-size: 22px; font-weight: 600; color: #111111;")
        header_layout.addWidget(title_label)
        
        header_layout.addStretch()
        
        # Center: Wavelength Controls (Pill shape)
        controls_frame = QFrame()
        controls_frame.setStyleSheet("""
            QFrame {
                background-color: #ffffff;
                border: 1px solid #e1e1e1;
                border-radius: 20px;
                padding: 4px;
            }
        """)
        controls_layout = QHBoxLayout(controls_frame)
        controls_layout.setContentsMargins(16, 4, 16, 4)
        controls_layout.setSpacing(10)
        
        wl_label = QLabel("Wavelength Range:")
        wl_label.setStyleSheet("font-size: 13px; font-weight: 500; color: #333333; border: none;")
        controls_layout.addWidget(wl_label)
        
        # Spinboxes
        spin_style = """
            QDoubleSpinBox {
                border: 1px solid #cccccc;
                border-radius: 4px;
                padding: 4px;
                min-width: 80px;
                background-color: #ffffff;
                color: #111111;
            }
        """
        self.wl_min_spin = QDoubleSpinBox()
        self.wl_min_spin.setRange(400.0, 900.0)
        self.wl_min_spin.setDecimals(1)
        self.wl_min_spin.setSingleStep(1.0)
        self.wl_min_spin.setValue(self.wavelength_min_nm)
        self.wl_min_spin.setSuffix(" nm")
        self.wl_min_spin.setStyleSheet(spin_style)
        
        self.wl_max_spin = QDoubleSpinBox()
        self.wl_max_spin.setRange(400.0, 900.0)
        self.wl_max_spin.setDecimals(1)
        self.wl_max_spin.setSingleStep(1.0)
        self.wl_max_spin.setValue(self.wavelength_max_nm)
        self.wl_max_spin.setSuffix(" nm")
        self.wl_max_spin.setStyleSheet(spin_style)
        
        controls_layout.addWidget(self.wl_min_spin)
        dash_label = QLabel("–")
        dash_label.setStyleSheet("color: #111111; font-weight: bold; border: none;")
        controls_layout.addWidget(dash_label)
        controls_layout.addWidget(self.wl_max_spin)
        
        self.apply_axis_btn = QPushButton("Apply Range")
        self.apply_axis_btn.setStyleSheet("""
            QPushButton {
                background-color: #f0f0f0;
                color: #111111;
                border: 1px solid #d0d0d0;
                border-radius: 4px;
                padding: 4px 12px;
                font-weight: 600;
                font-size: 12px;
            }
            QPushButton:hover { background-color: #e0e0e0; }
        """)
        self.apply_axis_btn.clicked.connect(self.on_wavelength_range_applied)
        controls_layout.addWidget(self.apply_axis_btn)
        
        header_layout.addWidget(controls_frame)
        header_layout.addStretch()
        
        # Right: Run Button
        self.run_exciton_btn = QPushButton("▶ Run Analysis")
        self.run_exciton_btn.setStyleSheet("""
            QPushButton {
                background-color: #111111;
                color: #ffffff;
                border: none;
                border-radius: 6px;
                padding: 10px 20px;
                font-weight: 600;
                font-size: 14px;
            }
            QPushButton:hover { background-color: #333333; }
            QPushButton:disabled { background-color: #cccccc; color: #888888; }
        """)
        self.run_exciton_btn.clicked.connect(self.run_exciton_distribution_calculation)
        self.run_exciton_btn.setEnabled(False)
        header_layout.addWidget(self.run_exciton_btn)
        
        main_layout.addWidget(header_container)
        
        # 2. Chart Grid (12 columns)
        grid_layout = QGridLayout()
        grid_layout.setSpacing(20)
        
        # Helper to create card frame
        def create_card():
            card = QFrame()
            card.setStyleSheet("""
                QFrame {
                    background-color: #ffffff;
                    border: 1px solid #e1e1e1;
                    border-radius: 12px;
                }
            """)
            return card

        # Left Card: Exciton Distribution (Span 7)
        left_card = create_card()
        left_layout = QVBoxLayout(left_card)
        left_layout.setContentsMargins(16, 16, 16, 16)
        
        # Card Header
        left_header = QHBoxLayout()
        left_title = QLabel("Exciton Distribution")
        left_title.setStyleSheet("font-size: 16px; font-weight: 600; border: none; color: #111111;")
        left_header.addWidget(left_title)
        left_header.addStretch()
        
        # Action buttons
        btn_style = """
            QPushButton { 
                border: 1px solid #e1e1e1; 
                border-radius: 4px; 
                padding: 4px 8px; 
                font-size: 11px; 
                background: white; 
                color: #333333;
            } 
            QPushButton:hover { background: #f9f9f9; }
        """
        details_btn = QPushButton("Details")
        details_btn.setStyleSheet(btn_style)
        export_btn = QPushButton("Export")
        export_btn.setStyleSheet(btn_style)
        left_header.addWidget(details_btn)
        left_header.addWidget(export_btn)
        
        left_layout.addLayout(left_header)
        
        # Content Area: Plot + Side Legend
        left_content = QHBoxLayout()
        left_content.setSpacing(0)
        
        # Plot Area (Stretch 3)
        self.exciton_distribution_plot = MiniPlot("Exciton Distribution")
        left_content.addWidget(self.exciton_distribution_plot, stretch=3)
        
        # Side Legend Area (Stretch 1)
        legend_container = QWidget()
        legend_container.setStyleSheet("background-color: #ffffff; border-left: 1px solid #f0f0f0;")
        legend_layout = QVBoxLayout(legend_container)
        legend_layout.setContentsMargins(0, 0, 0, 0)
        legend_layout.setSpacing(0)
        
        # Legend Header
        legend_header = QLabel("Pigments")
        legend_header.setStyleSheet("font-size: 12px; font-weight: 600; color: #666666; padding: 8px 12px; border-bottom: 1px solid #f0f0f0;")
        legend_layout.addWidget(legend_header)
        
        # Scrollable List
        self.pigment_scroll = QScrollArea()
        self.pigment_scroll.setWidgetResizable(True)
        self.pigment_scroll.setFrameShape(QFrame.NoFrame)
        self.pigment_scroll.setStyleSheet("""
            QScrollArea { background: transparent; border: none; }
            QScrollBar:vertical { width: 8px; background: #f0f0f0; }
            QScrollBar::handle:vertical { background: #d0d0d0; border-radius: 4px; }
        """)
        
        self.pigment_list_widget = QWidget()
        self.pigment_list_widget.setStyleSheet("background: transparent;")
        self.pigment_list_layout = QVBoxLayout(self.pigment_list_widget)
        self.pigment_list_layout.setContentsMargins(8, 8, 8, 8)
        self.pigment_list_layout.setSpacing(4)
        self.pigment_list_layout.addStretch()  # Push items to the top
        
        self.pigment_scroll.setWidget(self.pigment_list_widget)
        legend_layout.addWidget(self.pigment_scroll)
        
        left_content.addWidget(legend_container, stretch=1)
        left_layout.addLayout(left_content)
        
        grid_layout.addWidget(left_card, 0, 0, 2, 7) # Row 0, Col 0, RowSpan 2, ColSpan 7
        
        # Right Column (Span 5) - Single combined card
        right_card = create_card()
        right_layout = QVBoxLayout(right_card)
        right_layout.setContentsMargins(16, 16, 16, 16)
        
        right_header = QHBoxLayout()
        right_title = QLabel("Absorption & Exciton Analysis")
        right_title.setStyleSheet("font-size: 16px; font-weight: 600; border: none; color: #111111;")
        right_header.addWidget(right_title)
        right_header.addStretch()
        
        right_details = QPushButton("Details")
        right_details.setStyleSheet(btn_style)
        right_header.addWidget(right_details)
        
        right_layout.addLayout(right_header)
        
        self.combined_analysis_plot = MiniPlot("Combined Analysis")
        right_layout.addWidget(self.combined_analysis_plot)
        
        grid_layout.addWidget(right_card, 0, 7, 2, 5) # Row 0, Col 7, RowSpan 2, ColSpan 5
        
        # Set row stretches to be equal
        grid_layout.setRowStretch(0, 1)
        grid_layout.setRowStretch(1, 1)
        
        main_layout.addLayout(grid_layout)
        
        return tab_widget

    def connect_signals(self):
        """Connect all component signals"""

        # Tools Panel signals
        self.tools_panel.file_dropped.connect(self.on_file_dropped)
        self.tools_panel.open_project.connect(self.on_open_project)
        self.tools_panel.show_properties.connect(self.on_show_properties)
        self.tools_panel.parameters_changed.connect(self.on_parameters_changed)
        self.tools_panel.run_calculation.connect(self.on_run_calculation)
        self.tools_panel.view_results.connect(self.on_view_results)

        # Hamiltonian widget signals (Tab 1)
        self.hamiltonian_widget.run_hamiltonian_requested.connect(self.on_run_hamiltonian_requested)
        self.hamiltonian_widget.export_requested.connect(self.on_export_hamiltonian)
        self.hamiltonian_widget.view_full_matrix.connect(self.on_view_full_matrix)
        self.hamiltonian_widget.site_energy_updated.connect(self.on_site_energy_updated)
        self.hamiltonian_widget.parameters_changed.connect(self.on_hamiltonian_parameters_changed)
        self.hamiltonian_widget.domains_updated.connect(self.protein_viewer.update_domain_visualization)

        # Spectrum widget signals (Tab 2)
        self.spectrum_widget.run_spectra_requested.connect(self.on_run_spectra_requested)
        self.spectrum_widget.export_requested.connect(self.on_export_spectrum)
        self.spectrum_widget.parameters_changed.connect(self.on_spectra_parameters_changed)

        # Data table signals (Tab 3)
        # Data table removed from Data Analysis tab
        # self.data_table.row_selected.connect(self.on_pigment_selected)
        # self.data_table.export_requested.connect(self.on_export_table)
        # self.data_table.plot_requested.connect(self.on_plot_from_table)

    def on_file_dropped(self, file_path: str):
        """Handle file drop from drag-and-drop area"""
        if file_path:  # Non-empty path means actual file
            self.load_pdb_file(file_path)
        else:  # Empty path means user clicked, open dialog
            self.on_open_project()

    def on_export_hamiltonian(self):
        """Export Hamiltonian matrix to CSV"""
        if self.hamiltonian is None:
            QMessageBox.warning(self, "No Data", "No Hamiltonian available to export")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Hamiltonian Matrix",
            "hamiltonian.csv",
            "CSV Files (*.csv);;All Files (*)"
        )

        if file_path:
            try:
                np.savetxt(file_path, self.hamiltonian, delimiter=',', fmt='%.6f')
                self.status_message.emit(f"Hamiltonian exported to {Path(file_path).name}")
            except Exception as e:
                QMessageBox.critical(
                    self,
                    "Export Error",
                    f"Failed to export Hamiltonian:\n{str(e)}"
                )

    def on_view_full_matrix(self):
        """View full Hamiltonian matrix in dialog"""
        if self.hamiltonian is None:
            QMessageBox.warning(self, "No Data", "No Hamiltonian available")
            return

        from PyQt5.QtWidgets import QDialog, QVBoxLayout, QTextEdit, QDialogButtonBox

        dialog = QDialog(self)
        dialog.setWindowTitle("Full Hamiltonian Matrix")
        dialog.resize(700, 600)

        layout = QVBoxLayout(dialog)

        text_edit = QTextEdit()
        text_edit.setReadOnly(True)
        text_edit.setStyleSheet("background-color: #ffffff; color: #111111; font-family: monospace; font-size: 10px;")

        # Format matrix as text
        n = self.hamiltonian.shape[0]
        text = f"Hamiltonian Matrix ({n}×{n})\n"
        text += "=" * 80 + "\n\n"

        # Show header row with indices
        text += "     "
        for j in range(min(n, 10)):  # Show first 10 columns
            text += f"{j:>10d} "
        if n > 10:
            text += "..."
        text += "\n"
        text += "-" * 80 + "\n"

        # Show data rows
        for i in range(min(n, 20)):  # Show first 20 rows
            text += f"{i:3d}: "
            for j in range(min(n, 10)):
                text += f"{self.hamiltonian[i, j]:>10.2f} "
            if n > 10:
                text += "..."
            text += "\n"

        if n > 20:
            text += "...\n"

        text_edit.setPlainText(text)
        layout.addWidget(text_edit)

        # Buttons
        buttons = QDialogButtonBox(QDialogButtonBox.Ok)
        buttons.accepted.connect(dialog.accept)
        layout.addWidget(buttons)

        dialog.exec_()

    def on_export_spectrum(self):
        """Export spectrum data to CSV"""
        if self.spectrum_data is None:
            QMessageBox.warning(self, "No Data", "No spectrum available to export")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Spectrum",
            "spectrum.csv",
            "CSV Files (*.csv);;All Files (*)"
        )

        if file_path:
            try:
                import csv
                with open(file_path, 'w', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow(['Wavelength (nm)', 'Intensity'])
                    for wl, intensity in zip(self.spectrum_data['wavelengths'],
                                            self.spectrum_data['spectrum']):
                        writer.writerow([wl, intensity])
                self.status_message.emit(f"Spectrum exported to {Path(file_path).name}")
            except Exception as e:
                QMessageBox.critical(
                    self,
                    "Export Error",
                    f"Failed to export spectrum:\n{str(e)}"
                )

    def on_open_project(self):
        """Open and load a PDB file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Open PDB File",
            "",
            "PDB Files (*.pdb);;All Files (*)"
        )

        if file_path:
            self.load_pdb_file(file_path)

    def load_pdb_file(self, file_path: str):
        """Load and validate PDB file"""
        try:
            self.status_message.emit(f"Loading {Path(file_path).name}...")
            self.current_file = Path(file_path)

            # Load protein structure
            protein_structure = ProteinStructure.from_file_extended(
                file_path,
                self.current_file.name
            )
            protein_structure.pdb_path = file_path

            # Create pigment system
            self.pigment_system = PigmentSystem(protein_structure)

            if not self.pigment_system.pigments:
                raise ValueError("No pigments found in the structure")

            # Load TrEsp parameters for coupling calculations
            self.status_message.emit("Loading TrEsp parameters for coupling calculations...")
            tresp_mapping = {
                'CLA': 'CLA_IPPC',
                'CHL': 'CHL_IPPC',
            }

            n_loaded = 0
            for name, pigment in self.pigment_system.pigments.items():
                resname = pigment.get_resname()
                if resname in tresp_mapping:
                    tresp_param_name = tresp_mapping[resname]
                    cdc_param_name = resname  # CLA or CHL
                    try:
                        # Load BOTH CDC and TrEsp parameters for this pigment
                        pigment._load_calculation_params(tresp_dict_name=tresp_param_name, cdc_dict_name=cdc_param_name)
                        n_loaded += 1
                    except Exception as e:
                        logger.warning(f"Could not load TrEsp params for {name}: {e}")

            logger.info(f"Loaded TrEsp parameters for {n_loaded}/{len(self.pigment_system.pigments)} pigments")

            # Validate system
            validation = self.pigment_system.validate_system()

            # Count total atoms
            total_atoms = sum(1 for _ in protein_structure.pdb.get_atoms())

            # Update 3D viewer
            self.protein_viewer.load_structure(protein_structure, self.pigment_system)

            # Update tools panel
            self.tools_panel.update_project_info(
                name=self.current_file.name,
                pigments=validation.get('total_pigments', 0),
                atoms=total_atoms,
                validated=True
            )

            # Get current parameters
            params = self.get_all_parameters()

            # Initialize calculators
            self.initialize_calculators(params)

            # Get vacuum energies
            self.vacuum_energies = {}
            e0a = params.get('e0a', 14900.0)
            e0b = params.get('e0b', 15674.0)
            for pig_id, pigment in self.pigment_system.pigments.items():
                res_name = pigment.get_resname()
                if res_name == 'CLA':
                    self.vacuum_energies[pig_id] = e0a
                elif res_name == 'CHL':
                    self.vacuum_energies[pig_id] = e0b
                else:
                    self.vacuum_energies[pig_id] = e0a  # Default

            self.status_message.emit(f"Loaded {self.current_file.name} successfully")

            # Update results status
            self.tools_panel.update_result_status("site_energies", "pending")
            self.tools_panel.update_result_status("hamiltonian", "pending")
            self.tools_panel.update_result_status("spectrum", "pending")
            self.site_energy_overrides = {}
            self.calculated_site_energies = {}
            self.site_energies = {}
            self.pigment_labels = None
            self.hamiltonian_widget.set_run_enabled(True)
            self.spectrum_widget.set_run_enabled(False)

        except Exception as e:
            QMessageBox.critical(
                self,
                "Load Error",
                f"Failed to load PDB file:\n{str(e)}"
            )
            self.status_message.emit(f"Error: {str(e)}")

    def initialize_calculators(self, params: Dict[str, Any]):
        """Initialize all calculators with parameters"""
        if not self.pigment_system:
            return

        e0a = params.get('e0a', 14900.0)
        e0b = params.get('e0b', 15674.0)

        # Site energy calculator
        self.site_energy_calculator = SiteEnergyCalculator(
            dielectric_constant=params.get('dielectric_cdc', 2.0),
            e0a=e0a,
            e0b=e0b
        )

        # Hamiltonian calculator
        self.hamiltonian_calculator = HamiltonianCalculator(
            self.pigment_system,
            dielectric_cdc=params.get('dielectric_cdc', 2.0),
            dielectric_tresp=params.get('dielectric_tresp', 1.0),
            f_val=params.get('oscillator_strength', 0.72),
            E_0a=e0a,
            E_0b=e0b
        )

        # Pass hamiltonian calculator to widget for domain building
        self.hamiltonian_widget.set_hamiltonian_calculator(self.hamiltonian_calculator)
        self.hamiltonian_widget.E0a = e0a

        # Set coupling method
        self.hamiltonian_calculator.use_dipole_coupling = params.get('coupling_method') == 'dipole'

        # Spectrum calculator (will be created after diagonalization)
        self.spectrum_calculator = None

    def get_all_parameters(self) -> Dict[str, Any]:
        """Combine generic and tab-specific parameters"""
        params = self.tools_panel.get_parameters()
        params.update(self.hamiltonian_widget.get_parameters())
        params.update(self.spectrum_widget.get_parameters())
        return params

    def on_parameters_changed(self, params: Dict[str, Any]):
        """Handle parameter changes"""
        if self.pigment_system:
            combined = self.get_all_parameters()
            combined.update(params)
            self.initialize_calculators(combined)
            self.status_message.emit("Parameters updated")

    def on_hamiltonian_parameters_changed(self, params: Dict[str, Any]):
        """Handle Hamiltonian-specific parameter changes"""
        if self.hamiltonian_calculator:
            self.hamiltonian_calculator.use_dipole_coupling = params.get('coupling_method') == 'dipole'
        self.status_message.emit("Hamiltonian settings updated")

    def on_spectra_parameters_changed(self, params: Dict[str, Any]):
        """Handle spectra-specific parameter changes"""
        self.status_message.emit("Spectra settings updated")

    def on_run_hamiltonian_requested(self):
        """Handle Run Hamiltonian request from the Hamiltonian tab"""
        self.on_run_calculation("hamiltonian")

    def on_run_spectra_requested(self):
        """Handle Run Spectra request from the Spectra tab"""
        self.on_run_calculation("spectrum")

    def apply_site_energy_overrides(self, site_energies: Dict[str, float]) -> Dict[str, float]:
        """Apply user overrides to calculated site energies"""
        merged = dict(site_energies)
        for pig_id, energy in self.site_energy_overrides.items():
            if pig_id in merged:
                merged[pig_id] = energy
        return merged

    def extract_site_energies_from_hamiltonian(self) -> Dict[str, float]:
        """Extract site energies from the Hamiltonian diagonal"""
        if self.hamiltonian is None:
            return {}

        if hasattr(self.hamiltonian, "values"):
            diag = np.diag(self.hamiltonian.values)
            labels = list(self.hamiltonian.index)
        else:
            diag = np.diag(self.hamiltonian)
            labels = self.pigment_labels or []

        return {label: float(diag[i]) for i, label in enumerate(labels) if i < len(diag)}

    def apply_site_energies_to_hamiltonian(self, site_energies: Dict[str, float]):
        """Update Hamiltonian diagonal with the provided site energies"""
        if self.hamiltonian is None:
            return

        if hasattr(self.hamiltonian, "loc"):
            for pig_id, energy in site_energies.items():
                if pig_id in self.hamiltonian.index:
                    self.hamiltonian.loc[pig_id, pig_id] = energy
        else:
            if not self.pigment_labels:
                return
            for idx, pig_id in enumerate(self.pigment_labels):
                if pig_id in site_energies:
                    self.hamiltonian[idx, idx] = site_energies[pig_id]

    def refresh_hamiltonian_displays(self, update_interactions: bool = True):
        """Refresh Hamiltonian-dependent views after updates"""
        if self.hamiltonian is None or not self.pigment_labels:
            return

        ham_array = self.hamiltonian.values if hasattr(self.hamiltonian, "values") else self.hamiltonian
        self.eigenvalues, self.eigenvectors = np.linalg.eigh(ham_array)

        self.hamiltonian_widget.update_hamiltonian(self.hamiltonian, self.pigment_labels)
        self.hamiltonian_widget.update_eigenvalues(self.eigenvalues, self.eigenvectors)

        # Eigenvalue and interaction plots removed - replaced with exciton distribution plots

        if update_interactions:
            pass  # self.summary_cards.update_interactions(self.hamiltonian)

    def on_site_energy_updated(self, pigment_id: str, energy: float):
        """Handle user edits to site energy values"""
        self.site_energies[pigment_id] = energy
        baseline = self.calculated_site_energies.get(pigment_id)
        if baseline is not None and abs(energy - baseline) < 1e-6:
            self.site_energy_overrides.pop(pigment_id, None)
        else:
            self.site_energy_overrides[pigment_id] = energy

        if self.vacuum_energies:
            # self.data_table.update_data(self.site_energies, self.vacuum_energies)
            pass
        # self.summary_cards.update_statistics(self.site_energies)
        if self.site_energies:
            mean_energy = np.mean(list(self.site_energies.values()))
            self.tools_panel.update_result_status(
                "site_energies",
                "complete",
                f"Mean: {mean_energy:.1f} cm⁻¹"
            )

        if self.hamiltonian is not None:
            self.apply_site_energies_to_hamiltonian(self.site_energies)
            self.refresh_hamiltonian_displays(update_interactions=False)

        self.tools_panel.update_result_status("spectrum", "pending")
        self.status_message.emit(f"Updated site energy for {pigment_id}")

    def on_run_calculation(self, calc_type: str):
        """Run a calculation in background"""
        if not self.pigment_system:
            QMessageBox.warning(
                self,
                "No Project",
                "Please load a PDB file first"
            )
            return

        # Prevent multiple simultaneous calculations
        if self.current_worker and self.current_worker.isRunning():
            QMessageBox.warning(
                self,
                "Calculation Running",
                "Please wait for the current calculation to finish"
            )
            return

        try:
            if calc_type == "site_energies":
                self.run_site_energy_calculation()
            elif calc_type == "hamiltonian":
                self.run_hamiltonian_calculation()
            elif calc_type == "spectrum":
                self.run_spectrum_calculation()
            elif calc_type == "all":
                # Run complete workflow: site energies -> hamiltonian -> spectrum
                self.run_all_flag = True
                self.run_site_energy_calculation()
            else:
                QMessageBox.warning(
                    self,
                    "Unknown Calculation",
                    f"Unknown calculation type: {calc_type}"
                )
        except Exception as e:
            QMessageBox.critical(
                self,
                "Calculation Error",
                f"Failed to start calculation:\n{str(e)}"
            )

    def run_site_energy_calculation(self):
        """Calculate site energies in background"""
        self.tools_panel.update_result_status("site_energies", "running")
        self.status_message.emit("Calculating site energies...")

        self.current_worker = CalculationWorker(
            "site_energies",
            self.site_energy_calculator,
            pigment_system=self.pigment_system
        )
        self.current_worker.finished.connect(self.on_calculation_finished)
        self.current_worker.error.connect(self.on_calculation_error)
        self.current_worker.progress.connect(self.on_calculation_progress)
        self.current_worker.start()

    def run_hamiltonian_calculation(self):
        """Calculate Hamiltonian in background"""
        self.tools_panel.update_result_status("hamiltonian", "running")
        self.status_message.emit("Constructing Hamiltonian...")
        self.spectrum_widget.set_run_enabled(False)

        coupling_method = self.hamiltonian_widget.get_parameters().get('coupling_method', 'tresp')
        self.current_worker = CalculationWorker(
            "hamiltonian",
            self.hamiltonian_calculator,
            coupling_calc=coupling_method
        )
        self.current_worker.finished.connect(self.on_calculation_finished)
        self.current_worker.error.connect(self.on_calculation_error)
        self.current_worker.progress.connect(self.on_calculation_progress)
        self.current_worker.start()

    def run_spectrum_calculation(self):
        """Calculate spectrum using advanced Renger method"""
        if self.hamiltonian is None or self.eigenvalues is None:
            QMessageBox.warning(
                self,
                "Missing Data",
                "Please calculate Hamiltonian first"
            )
            return

        if self.site_energies:
            self.apply_site_energies_to_hamiltonian(self.site_energies)
            self.refresh_hamiltonian_displays(update_interactions=False)

        # Get all parameters including spectrum-specific ones (n_ensemble, temperature, etc.)
        params = self.get_all_parameters()

        # Get domains from Hamiltonian widget
        domains = self.hamiltonian_widget.get_current_domains()

        # Create spectra calculator with Renger lineshape theory
        from Alprotein.calculators.spectra_calculator import SpectraCalculator
        spectra_calculator = SpectraCalculator(
            temperature=params.get('temperature', 77.0),
            disorder_sigma=120.0,  # FWHM in cm⁻¹
            n_ensemble=params.get('n_ensemble', 100),
            dt=0.1,  # fs
            t_max=2000.0  # fs
        )

        # Store spectra_calculator for exciton distribution
        self.spectra_calculator = spectra_calculator

        self.tools_panel.update_result_status("spectrum", "running")
        self.status_message.emit("Calculating absorption spectrum (Renger lineshape theory)...")

        # Run advanced Renger calculation in background
        self.current_worker = CalculationWorker(
            "spectrum_renger",
            spectra_calculator,
            hamiltonian=self.hamiltonian,
            domains=domains,
            pigment_system=self.pigment_system,
            include_vibronic=params.get('include_vibronic', True),
            use_ensemble=params.get('use_ensemble', True)
        )
        self.current_worker.finished.connect(self.on_calculation_finished)
        self.current_worker.error.connect(self.on_calculation_error)
        self.current_worker.progress.connect(self.on_calculation_progress)
        self.current_worker.start()

    def run_exciton_distribution_calculation(self):
        """Calculate exciton distributions"""
        if self.hamiltonian is None:
            QMessageBox.warning(
                self,
                "Missing Data",
                "Please calculate Hamiltonian first"
            )
            return

        # Get domains from Hamiltonian widget
        domains = self.hamiltonian_widget.get_current_domains()

        # Create exciton calculator (share parameters with spectra calculator)
        if hasattr(self, 'spectra_calculator') and self.spectra_calculator:
            self.exciton_calculator = ExcitonCalculator(
                disorder_sigma=self.spectra_calculator.disorder_sigma,
                n_ensemble=self.spectra_calculator.n_ensemble,
                temperature=self.spectra_calculator.temperature
            )
        else:
            # Fallback if spectra not calculated yet
            params = self.get_all_parameters()
            self.exciton_calculator = ExcitonCalculator(
                disorder_sigma=55.0,  # Default sigma (120 FWHM converted)
                n_ensemble=params.get('n_ensemble', 100),
                temperature=params.get('temperature', 77.0)
            )

        # Update tools panel status
        self.tools_panel.update_result_status("exciton_distribution", "running")

        self.status_message.emit("Calculating exciton distributions...")

        # Run exciton distribution calculation in background
        self.current_worker = CalculationWorker(
            "exciton_distribution",
            self.exciton_calculator,
            hamiltonian=self.hamiltonian,
            domains=domains,
            pigment_system=self.pigment_system
        )
        self.current_worker.finished.connect(self.on_calculation_finished)
        self.current_worker.error.connect(self.on_calculation_error)
        self.current_worker.progress.connect(self.on_calculation_progress)
        self.current_worker.start()

    def on_calculation_finished(self, calc_type: str, result: Any):
        """Handle calculation completion and distribute data to all tabs"""
        try:
            logger.info(f"=== on_calculation_finished called with calc_type: '{calc_type}' ===")

            if calc_type == "site_energies":
                self.calculated_site_energies = result
                self.site_energies = self.apply_site_energy_overrides(result)

                # Update Tab 1: Hamiltonian widget - Site energies table
                pigment_labels = {pig_id: pig.get_resname()
                                 for pig_id, pig in self.pigment_system.pigments.items()}
                self.hamiltonian_widget.update_site_energies_table(
                    self.site_energies, pigment_labels
                )

                # Update Tab 3: Data table (removed)
                # self.data_table.update_data(
                #     self.site_energies,
                #     self.vacuum_energies
                # )

                # Update Tab 3: Summary cards - REMOVED
                # self.summary_cards.update_statistics(self.site_energies)
                # pig_counts = {}
                # for pig in self.pigment_system.pigments.values():
                #     res = pig.get_resname()
                #     pig_counts[res] = pig_counts.get(res, 0) + 1
                # self.summary_cards.update_pigments(pig_counts)

                # Calculate mean energy for status
                mean_energy = np.mean(list(self.site_energies.values()))
                self.tools_panel.update_result_status(
                    "site_energies",
                    "complete",
                    f"Mean: {mean_energy:.1f} cm⁻¹"
                )

                if self.hamiltonian is not None:
                    self.apply_site_energies_to_hamiltonian(self.site_energies)
                    self.refresh_hamiltonian_displays(update_interactions=False)

                self.status_message.emit("Site energies calculated successfully")

            elif calc_type == "hamiltonian":
                self.hamiltonian = result
                if hasattr(self.hamiltonian, "index"):
                    self.pigment_labels = list(self.hamiltonian.index)
                else:
                    self.pigment_labels = list(self.pigment_system.pigments.keys())

                self.calculated_site_energies = self.extract_site_energies_from_hamiltonian()
                self.site_energies = self.apply_site_energy_overrides(self.calculated_site_energies)
                self.apply_site_energies_to_hamiltonian(self.site_energies)

                pigment_labels = {pig_id: pig.get_resname()
                                 for pig_id, pig in self.pigment_system.pigments.items()}
                self.hamiltonian_widget.update_site_energies_table(
                    self.site_energies, pigment_labels
                )

                # self.data_table.update_data(
                #     self.site_energies,
                #     self.vacuum_energies
                # )
                # self.summary_cards.update_statistics(self.site_energies)

                self.refresh_hamiltonian_displays(update_interactions=True)

                # Get coupling range
                # Convert to numpy array if it's a DataFrame
                ham_array = self.hamiltonian.values if hasattr(self.hamiltonian, 'values') else self.hamiltonian
                n = ham_array.shape[0]
                off_diag = []
                for i in range(n):
                    for j in range(i+1, n):
                        off_diag.append(abs(ham_array[i, j]))

                if off_diag:
                    coupling_range = f"{min(off_diag):.1f} - {max(off_diag):.1f} cm⁻¹"
                else:
                    coupling_range = "N/A"

                if self.site_energies:
                    mean_energy = np.mean(list(self.site_energies.values()))
                    self.tools_panel.update_result_status(
                        "site_energies",
                        "complete",
                        f"Mean: {mean_energy:.1f} cm⁻¹"
                    )

                self.tools_panel.update_result_status(
                    "hamiltonian",
                    "complete",
                    coupling_range
                )
                self.spectrum_widget.set_run_enabled(True)

                self.status_message.emit("Hamiltonian constructed and diagonalized")

            elif calc_type == "spectrum_renger" or calc_type == "spectrum":
                self.spectrum_data = result
                method = result.get('method', 'unknown')

                # Check if we have the new format (with fluorescence) or old format
                if 'wavelengths_abs' in result and 'fluorescence' in result:
                    # New format: Update both absorption and fluorescence
                    wavelengths_abs = result['wavelengths_abs']
                    absorption = result['absorption']
                    wavelengths_fl = result['wavelengths_fl']
                    fluorescence = result['fluorescence']

                    self.spectrum_widget.update_spectra(
                        wavelengths_abs=wavelengths_abs,
                        absorption=absorption,
                        wavelengths_fl=wavelengths_fl,
                        fluorescence=fluorescence
                    )

                    # Find absorption peak for status
                    peak_idx = np.argmax(absorption)
                    peak_wl = wavelengths_abs[peak_idx]

                    # Find fluorescence peak
                    fl_peak_idx = np.argmax(fluorescence)
                    fl_peak_wl = wavelengths_fl[fl_peak_idx]

                    method_label = "Renger" if method == "renger" else "Simple"
                    self.tools_panel.update_result_status(
                        "spectrum",
                        "complete",
                        f"Abs: {peak_wl:.1f} nm, Fl: {fl_peak_wl:.1f} nm ({method_label})"
                    )

                    self.status_message.emit(
                        f"Absorption and fluorescence spectra calculated successfully ({method_label} method)"
                    )

                    # Enable exciton distribution button
                    logger.info("Enabling exciton distribution button (new format)")
                    self.run_exciton_btn.setEnabled(True)
                    self.run_exciton_btn.update()  # Force GUI refresh
                    logger.info(f"Button enabled state: {self.run_exciton_btn.isEnabled()}")

                else:
                    # Old format: Only absorption (for backward compatibility)
                    wavelengths = result['wavelengths']
                    spectrum = result['spectrum']
                    A_00 = result.get('A_00', None)
                    A_01 = result.get('A_01', None)

                    self.spectrum_widget.update_spectrum_with_components(
                        wavelengths=wavelengths,
                        spectrum=spectrum,
                        A_00=A_00,
                        A_01=A_01
                    )

                    peak_idx = np.argmax(spectrum)
                    peak_wl = wavelengths[peak_idx]

                    method_label = "Renger" if method == "renger" else "Simple"
                    self.tools_panel.update_result_status(
                        "spectrum",
                        "complete",
                        f"Peak: {peak_wl:.1f} nm ({method_label})"
                    )

                    self.status_message.emit(f"Spectrum calculated successfully ({method_label} method)")

                    # Enable exciton distribution button
                    logger.info("Enabling exciton distribution button")
                    self.run_exciton_btn.setEnabled(True)
                    self.run_exciton_btn.update()  # Force GUI refresh
                    logger.info(f"Button enabled state: {self.run_exciton_btn.isEnabled()}")

            elif calc_type == "exciton_distribution":
                # Store exciton distributions
                exciton_labels = result['labels']
                self.exciton_distributions = result['distributions']
                self.exciton_labels = exciton_labels

                # Update tools panel status
                self.tools_panel.update_result_status(
                    "exciton_distribution",
                    "complete",
                    f"{len(exciton_labels)} pigments"
                )

                self.status_message.emit("Exciton distributions calculated successfully")

                # Update plots in Data Analysis tab
                self.update_exciton_plots(self.exciton_distributions, exciton_labels)

            # Auto-chain calculations if running complete workflow
            if self.run_all_flag:
                if calc_type == "site_energies":
                    # After site energies, run Hamiltonian
                    self.status_message.emit("Auto-starting Hamiltonian calculation...")
                    self.run_hamiltonian_calculation()
                elif calc_type == "hamiltonian":
                    # After Hamiltonian, run spectrum
                    self.status_message.emit("Auto-starting spectrum calculation...")
                    self.run_spectrum_calculation()
                elif calc_type == "spectrum_renger" or calc_type == "spectrum":
                    # After spectrum, run exciton distribution
                    self.status_message.emit("Auto-starting exciton distribution calculation...")
                    self.run_exciton_distribution_calculation()
                elif calc_type == "exciton_distribution":
                    # Complete workflow finished
                    self.run_all_flag = False
                    self.status_message.emit("Complete workflow finished successfully!")
                    QMessageBox.information(
                        self,
                        "Workflow Complete",
                        "All calculations completed successfully!\n\n"
                        "✓ Site Energies\n"
                        "✓ Hamiltonian\n"
                        "✓ Absorption & Fluorescence Spectra\n"
                        "✓ Exciton Distributions\n\n"
                        "View results in all tabs."
                    )

        except Exception as e:
            self.run_all_flag = False  # Reset flag on error
            self.on_calculation_error(calc_type, str(e))

    def on_calculation_error(self, calc_type: str, error_msg: str):
        """Handle calculation error"""
        # Reset workflow flag if active
        self.run_all_flag = False

        # Log the error with type information
        logger.error(f"Calculation error in {calc_type}:")
        logger.error(f"  Error message type: {type(error_msg)}")
        logger.error(f"  Error message: {error_msg}")

        self.tools_panel.update_result_status(calc_type, "error", error_msg[:50])
        self.status_message.emit(f"Error in {calc_type}: {error_msg}")

        QMessageBox.critical(
            self,
            "Calculation Error",
            f"Error in {calc_type}:\n{error_msg}"
        )

    def on_calculation_progress(self, message: str, percentage: int):
        """Handle calculation progress"""
        self.progress_update.emit(message, percentage)

    def on_show_properties(self):
        """Show project properties dialog"""
        if not self.pigment_system:
            return

        # Create simple properties message
        props = []
        props.append(f"File: {self.current_file.name if self.current_file else 'N/A'}")
        props.append(f"Total Pigments: {len(self.pigment_system.pigments)}")

        # Get atom count from BioPython structure
        atom_count = len(list(self.pigment_system.structure.pdb.get_atoms()))
        props.append(f"Total Atoms: {atom_count}")
        props.append("\nPigment Types:")

        # Count pigment types
        types = {}
        for pig_id, pig in self.pigment_system.pigments.items():
            res = pig.get_resname()
            types[res] = types.get(res, 0) + 1

        for res, count in types.items():
            props.append(f"  {res}: {count}")

        QMessageBox.information(
            self,
            "Project Properties",
            "\n".join(props)
        )

    def on_view_results(self, result_type: str):
        """View all results or specific result"""
        # For now, just show a message
        # In full implementation, this would open detailed views
        self.status_message.emit(f"Viewing {result_type} results")

    def on_view_details(self, plot_type: str):
        """View detailed plot (full size)"""
        # For now, just show a message
        # In full implementation, this would open full-size plot windows
        self.status_message.emit(f"Opening detailed {plot_type} view")

    def on_pigment_selected(self, pigment_id: str):
        """Handle pigment selection from table"""
        # Highlight in 3D viewer
        self.protein_viewer.highlight_pigment(pigment_id)
        self.status_message.emit(f"Selected pigment: {pigment_id}")

    def on_export_table(self):
        """Export data table"""
        if not self.site_energies:
            QMessageBox.warning(
                self,
                "No Data",
                "No data available to export"
            )
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Data",
            "",
            "CSV Files (*.csv);;JSON Files (*.json)"
        )

        if file_path:
            try:
                # The data table widget handles the actual export
                self.status_message.emit(f"Exported to {Path(file_path).name}")
            except Exception as e:
                QMessageBox.critical(
                    self,
                    "Export Error",
                    f"Failed to export:\n{str(e)}"
                )

    def on_plot_from_table(self):
        """Create plot from table data"""
        # For now, just update the dashboard
        if self.site_energies and self.vacuum_energies:
            self.analysis_dashboard.update_site_energies(
                self.site_energies,
                self.vacuum_energies
            )
            self.status_message.emit("Updated plots from data")

    def on_wavelength_range_applied(self):
        """Apply user-selected wavelength axis limits to exciton plots"""
        min_val = self.wl_min_spin.value()
        max_val = self.wl_max_spin.value()

        if min_val >= max_val:
            QMessageBox.warning(
                self,
                "Invalid Range",
                "Minimum wavelength must be smaller than maximum."
            )
            return

        self.wavelength_min_nm = min_val
        self.wavelength_max_nm = max_val

        if self.exciton_distributions and self.exciton_labels:
            self.update_exciton_plots(self.exciton_distributions, self.exciton_labels)
        else:
            self.status_message.emit(
                f"Wavelength range set to {min_val:.1f}–{max_val:.1f} nm (plots will update after calculation)"
            )

    def update_exciton_plots(self, distributions, labels):
        """Update exciton distribution plots in Data Analysis tab with Gaussian curves"""
        try:
            import matplotlib.pyplot as plt
            import numpy as np
            from scipy.interpolate import interp1d
            from scipy.ndimage import gaussian_filter1d

            # Persist labels for later range adjustments
            self.exciton_labels = labels

            x_min, x_max = self.wavelength_min_nm, self.wavelength_max_nm

            # Debug logging
            logger.info(f"update_exciton_plots called with {len(labels)} labels")

            # Plot 1: Exciton Distribution for all pigments (Left Large Card)
            fig1 = self.exciton_distribution_plot.figure
            fig1.clear()
            ax1 = fig1.add_subplot(111)

            # Plot 2: Combined Analysis (Right Card)
            fig2 = self.combined_analysis_plot.figure
            fig2.clear()
            
            # Create subplots sharing x-axis
            gs = fig2.add_gridspec(2, 1, height_ratios=[1, 1.2], hspace=0.05)
            ax_top = fig2.add_subplot(gs[0])
            ax_bottom = fig2.add_subplot(gs[1], sharex=ax_top)

            # Store line objects for interactive legend
            lines_left = []
            lines_right = []
            line_labels = []

            # Common settings
            num_pigments_to_show = 10
            colors = plt.cm.tab20(np.linspace(0, 1, len(labels)))

            # Plot each pigment's distribution as smooth Gaussian curve
            for i, label in enumerate(labels[:num_pigments_to_show]):
                energies, probs = distributions[label]
                energies = np.array(energies)
                probs = np.array(probs)
                
                # Convert energies to wavelengths
                wavelengths = np.array([1e7 / e if e > 0 else 0 for e in energies])
                
                # Remove any invalid values
                valid_mask = (wavelengths > 0) & (probs > 0) & np.isfinite(wavelengths) & np.isfinite(probs)
                wavelengths = wavelengths[valid_mask]
                probs = probs[valid_mask]
                
                if len(wavelengths) < 3:
                    # Add dummy lines to keep indices aligned if we skip
                    # But better to just continue and handle alignment carefully
                    # For simplicity in this fix, we'll just skip and not add to lists
                    # which might cause index mismatch if we relied on 'i'. 
                    # But we append to lists, so list indices will match.
                    continue
                
                # Create histogram bins for the wavelength range
                bins = np.linspace(x_min, x_max, 200)
                
                # Create weighted histogram (probability density)
                hist, bin_edges = np.histogram(wavelengths, bins=bins, weights=probs, density=False)
                bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
                
                # Normalize histogram
                if hist.sum() > 0:
                    hist = hist / hist.sum()
                
                # Apply Gaussian smoothing for smooth curves
                hist_smooth = gaussian_filter1d(hist, sigma=3)
                
                # Plot on Left (ax1)
                line1, = ax1.plot(bin_centers, hist_smooth, linewidth=2.5, 
                                label=label.split('_')[-1], color=colors[i], alpha=0.8)
                lines_left.append(line1)
                
                # Plot on Right Bottom (ax_bottom)
                line2, = ax_bottom.plot(bin_centers, hist_smooth, linewidth=2.0, 
                                     label=label.split('_')[-1], color=colors[i], alpha=0.9)
                lines_right.append(line2)
                
                line_labels.append(label.split('_')[-1])

            # Styling for Left Plot
            ax1.set_xlabel('Wavelength (nm)', fontsize=10, fontweight='bold')
            ax1.set_ylabel('Density', fontsize=10, fontweight='bold')
            ax1.tick_params(axis='both', which='major', labelsize=9)
            ax1.set_xlim(x_min, x_max)
            ax1.grid(True, alpha=0.3, linestyle='--')
            fig1.subplots_adjust(bottom=0.15, top=0.95, left=0.15, right=0.98)
            
            # Styling for Right Top (Absorption)
            if self.spectrum_data is not None:
                wavelengths_abs = self.spectrum_data.get('wavelengths_abs',
                                                        self.spectrum_data.get('wavelengths'))
                absorption = self.spectrum_data.get('absorption',
                                                   self.spectrum_data.get('spectrum'))
                ax_top.plot(wavelengths_abs, absorption, 'k-', linewidth=3.0, label='Absorption')
                ax_top.set_ylabel('Absorption (a.u.)', fontsize=10, fontweight='bold')
                ax_top.tick_params(axis='both', which='major', labelsize=9)
                ax_top.grid(False)
                plt.setp(ax_top.get_xticklabels(), visible=False)
                ax_top.spines['top'].set_visible(False)
                ax_top.spines['right'].set_visible(False)
            else:
                ax_top.text(0.5, 0.5, "No Spectrum Data", ha='center', va='center')
                ax_top.axis('off')

            # Styling for Right Bottom (Exciton)
            ax_bottom.set_xlabel('Wavelength (nm)', fontsize=10, fontweight='bold')
            ax_bottom.set_ylabel('Density', fontsize=10, fontweight='bold')
            ax_bottom.tick_params(axis='both', which='major', labelsize=9)
            ax_bottom.set_xlim(x_min, x_max)
            ax_bottom.spines['top'].set_visible(False)
            ax_bottom.spines['right'].set_visible(False)
            ax_bottom.grid(axis='x', alpha=0.3)
            fig2.subplots_adjust(bottom=0.15, top=0.95, left=0.18, right=0.95)

            # Populate the side pigment list (Controls BOTH plots)
            while self.pigment_list_layout.count() > 1: 
                item = self.pigment_list_layout.takeAt(0)
                if item.widget():
                    item.widget().deleteLater()
            
            for i, (line1, line2, label) in enumerate(zip(lines_left, lines_right, line_labels)):
                # Get color from line
                color_rgba = line1.get_color()
                if isinstance(color_rgba, (tuple, list, np.ndarray)):
                    r, g, b = [int(c * 255) for c in color_rgba[:3]]
                    hex_color = f"#{r:02x}{g:02x}{b:02x}"
                else:
                    hex_color = color_rgba
                
                item_widget = QWidget()
                item_layout = QHBoxLayout(item_widget)
                item_layout.setContentsMargins(4, 2, 4, 2)
                item_layout.setSpacing(8)
                
                color_box = QLabel()
                color_box.setFixedSize(12, 12)
                color_box.setStyleSheet(f"background-color: {hex_color}; border-radius: 6px;")
                item_layout.addWidget(color_box)
                
                chk = QCheckBox(label)
                chk.setChecked(True)
                chk.setStyleSheet("font-size: 11px; color: #333333;")
                
                # Connect checkbox to BOTH lines
                def toggle_lines(state, l1=line1, l2=line2):
                    visible = (state == Qt.Checked)
                    l1.set_visible(visible)
                    l2.set_visible(visible)
                    fig1.canvas.draw()
                    fig2.canvas.draw()
                    
                chk.stateChanged.connect(toggle_lines)
                item_layout.addWidget(chk)
                item_layout.addStretch()
                
                self.pigment_list_layout.insertWidget(i, item_widget)

            # Store lines for interactive toggling (via plot click)
            self.exciton_distribution_plot.legend_lines = lines_left
            self.exciton_distribution_plot.legend_labels = line_labels
            
            # Make lines interactive
            for line in lines_left:
                line.set_picker(5)
            
            # Connect click event on Left Plot
            def on_plot_pick(event):
                # We need to handle the pick manually to update both plots
                if event.artist in lines_left:
                    idx = lines_left.index(event.artist)
                    line1 = lines_left[idx]
                    line2 = lines_right[idx]
                    
                    # Toggle
                    visible = not line1.get_visible()
                    line1.set_visible(visible)
                    line2.set_visible(visible)
                    
                    # Update styles
                    if visible:
                        line1.set_alpha(0.8); line1.set_linewidth(2.5)
                        line2.set_alpha(0.9); line2.set_linewidth(2.0)
                    else:
                        line1.set_alpha(0.2); line1.set_linewidth(1.0)
                        line2.set_alpha(0.2); line2.set_linewidth(1.0)
                    
                    fig1.canvas.draw()
                    fig2.canvas.draw()
                    
                    # Sync Checkbox
                    item = self.pigment_list_layout.itemAt(idx)
                    if item and item.widget():
                        chk_widget = item.widget().layout().itemAt(1).widget()
                        if isinstance(chk_widget, QCheckBox):
                            chk_widget.blockSignals(True)
                            chk_widget.setChecked(visible)
                            chk_widget.blockSignals(False)

            fig1.canvas.mpl_connect('pick_event', on_plot_pick)

            self.exciton_distribution_plot.canvas.draw()
            self.exciton_distribution_plot.update()
            self.combined_analysis_plot.canvas.draw()
            self.combined_analysis_plot.update()

            logger.info("Exciton distribution plots updated successfully")
            self.status_message.emit("Exciton distribution plots updated")

        except Exception as e:
            import traceback
            error_details = traceback.format_exc()
            logger.error(f"Failed to update exciton plots: {e}\n{error_details}")
            self.status_message.emit(f"Error updating exciton plots: {str(e)}")

    def on_line_pick(self, event, figure, lines, line_labels):
        """Handle line click to toggle line visibility"""
        try:
            # Get the line that was clicked
            if event.artist in lines:
                line = event.artist
                idx = lines.index(line)
                
                # Toggle visibility
                visible = not line.get_visible()
                line.set_visible(visible)
                
                # Change line appearance
                if visible:
                    line.set_alpha(0.8)
                    line.set_linewidth(2.5)
                else:
                    line.set_alpha(0.2)
                    line.set_linewidth(1.0)
                
                # Redraw
                figure.canvas.draw()
                figure.canvas.flush_events()
                
                label = line_labels[idx] if idx < len(line_labels) else "line"
                status = "shown" if visible else "hidden"
                self.status_message.emit(f"Peak {label} {status}")
                    
        except Exception as e:
            logger.error(f"Error in line pick handler: {e}")

    def on_legend_pick(self, event, plot_widget):
        """Handle legend click to toggle line visibility (deprecated - using on_line_pick instead)"""
        try:
            # Get the legend line that was clicked
            legline = event.artist
            
            # Find which line it corresponds to
            if hasattr(plot_widget, 'legend') and hasattr(plot_widget, 'lines'):
                legend = plot_widget.legend
                legend_lines = legend.get_lines()
                
                # Find index of clicked legend line
                try:
                    idx = list(legend_lines).index(legline)
                    origline = plot_widget.lines[idx]
                    
                    # Toggle visibility
                    visible = not origline.get_visible()
                    origline.set_visible(visible)
                    
                    # Change legend line appearance
                    if visible:
                        legline.set_alpha(1.0)
                        legline.set_linewidth(2.0)
                    else:
                        legline.set_alpha(0.2)
                        legline.set_linewidth(1.0)
                    
                    # Redraw
                    plot_widget.canvas.draw()
                    plot_widget.canvas.flush_events()
                    
                    label = plot_widget.line_labels[idx] if idx < len(plot_widget.line_labels) else "line"
                    status = "shown" if visible else "hidden"
                    self.status_message.emit(f"Peak {label} {status}")
                    
                except (ValueError, IndexError) as e:
                    logger.warning(f"Could not find legend line: {e}")
                    
        except Exception as e:
            logger.error(f"Error in legend pick handler: {e}")
            import traceback
            traceback.print_exc()

    def clear_all_data(self):
        """Clear all loaded data and results"""
        self.pigment_system = None
        self.site_energies = {}
        self.calculated_site_energies = {}
        self.site_energy_overrides = {}
        self.vacuum_energies = {}
        self.hamiltonian = None
        self.eigenvalues = None
        self.eigenvectors = None
        self.spectrum_data = None
        self.exciton_distributions = None
        self.exciton_labels = None
        self.pigment_labels = None

        # self.data_table.clear_data()  # Data table removed
        self.protein_viewer.clear_structure()
        self.hamiltonian_widget.clear()
        self.spectrum_widget.clear()
        self.hamiltonian_widget.set_run_enabled(False)
        self.spectrum_widget.set_run_enabled(False)

        # Clear exciton plots
        self.exciton_distribution_plot.clear_plot()
        self.combined_analysis_plot.clear_plot()

        # Disable exciton distribution button
        self.run_exciton_btn.setEnabled(False)

        self.status_message.emit("All data cleared")
