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
    QMessageBox, QFileDialog, QTabWidget
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
from Alprotein.gui.widgets.summary_cards_widget import SummaryCardsWidget
from Alprotein.core.pigment_system import PigmentSystem
from Alprotein.core.protein_structure import ProteinStructure
from Alprotein.calculators.site_energy_calculator import SiteEnergyCalculator
from Alprotein.calculators.hamiltonian_calculator import HamiltonianCalculator


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

        # Results storage
        self.site_energies: Dict[str, float] = {}
        self.calculated_site_energies: Dict[str, float] = {}
        self.site_energy_overrides: Dict[str, float] = {}
        self.vacuum_energies: Dict[str, float] = {}
        self.hamiltonian: Optional[np.ndarray] = None
        self.eigenvalues: Optional[np.ndarray] = None
        self.eigenvectors: Optional[np.ndarray] = None
        self.spectrum_data: Optional[Dict[str, np.ndarray]] = None
        self.pigment_labels: Optional[List[str]] = None

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
        self.workspace_tabs.setStyleSheet("""
            QTabWidget::pane {
                border: 1px solid #e1e1e1;
                border-radius: 2px;
                background-color: #ffffff;
            }
            QTabBar {
                qproperty-expanding: 0;
            }
            QTabBar::tab {
                background-color: #f3f3f3;
                color: #111111;
                padding: 10px 20px;
                min-width: 120px;
                border: 1px solid #e1e1e1;
                border-bottom: none;
                border-top-left-radius: 2px;
                border-top-right-radius: 2px;
                margin-right: 2px;
                font-size: 13px;
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
        """Create enhanced Data Analysis tab with summary cards"""
        tab_widget = QWidget()
        tab_layout = QVBoxLayout(tab_widget)
        tab_layout.setContentsMargins(16, 16, 16, 16)
        tab_layout.setSpacing(16)

        # Summary cards at top
        self.summary_cards = SummaryCardsWidget()
        tab_layout.addWidget(self.summary_cards, stretch=0)

        # Data table in middle
        self.data_table = DataTableWidget()
        tab_layout.addWidget(self.data_table, stretch=2)

        # Bottom plots (2 side-by-side)
        plots_widget = QWidget()
        plots_layout = QHBoxLayout(plots_widget)
        plots_layout.setSpacing(16)

        # Left plot: Eigenvalue distribution
        self.eigenvalue_plot = MiniPlot("Eigenvalue Distribution")
        plots_layout.addWidget(self.eigenvalue_plot)

        # Right plot: Interaction plot
        self.interaction_plot = MiniPlot("Coupling Interactions")
        plots_layout.addWidget(self.interaction_plot)

        tab_layout.addWidget(plots_widget, stretch=1)

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
        self.data_table.row_selected.connect(self.on_pigment_selected)
        self.data_table.export_requested.connect(self.on_export_table)
        self.data_table.plot_requested.connect(self.on_plot_from_table)

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

        self.eigenvalue_plot.plot_histogram(
            self.eigenvalues,
            "Energy (cm⁻¹)",
            "Eigenvalue Distribution"
        )

        if update_interactions:
            self.summary_cards.update_interactions(self.hamiltonian)
            self.interaction_plot.plot_heatmap(self.hamiltonian, "Coupling Matrix")

    def on_site_energy_updated(self, pigment_id: str, energy: float):
        """Handle user edits to site energy values"""
        self.site_energies[pigment_id] = energy
        baseline = self.calculated_site_energies.get(pigment_id)
        if baseline is not None and abs(energy - baseline) < 1e-6:
            self.site_energy_overrides.pop(pigment_id, None)
        else:
            self.site_energy_overrides[pigment_id] = energy

        if self.vacuum_energies:
            self.data_table.update_data(self.site_energies, self.vacuum_energies)
        self.summary_cards.update_statistics(self.site_energies)
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

    def on_calculation_finished(self, calc_type: str, result: Any):
        """Handle calculation completion and distribute data to all tabs"""
        try:
            if calc_type == "site_energies":
                self.calculated_site_energies = result
                self.site_energies = self.apply_site_energy_overrides(result)

                # Update Tab 1: Hamiltonian widget - Site energies table
                pigment_labels = {pig_id: pig.get_resname()
                                 for pig_id, pig in self.pigment_system.pigments.items()}
                self.hamiltonian_widget.update_site_energies_table(
                    self.site_energies, pigment_labels
                )

                # Update Tab 3: Data table
                self.data_table.update_data(
                    self.site_energies,
                    self.vacuum_energies
                )

                # Update Tab 3: Summary cards - Statistical Summary
                self.summary_cards.update_statistics(self.site_energies)

                # Update Tab 3: Summary cards - Pigment Counts
                pig_counts = {}
                for pig in self.pigment_system.pigments.values():
                    res = pig.get_resname()
                    pig_counts[res] = pig_counts.get(res, 0) + 1
                self.summary_cards.update_pigments(pig_counts)

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

                self.data_table.update_data(
                    self.site_energies,
                    self.vacuum_energies
                )
                self.summary_cards.update_statistics(self.site_energies)

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
                    # Complete workflow finished
                    self.run_all_flag = False
                    self.status_message.emit("Complete workflow finished successfully!")
                    QMessageBox.information(
                        self,
                        "Workflow Complete",
                        "All calculations completed successfully!\n\n"
                        "✓ Site Energies\n"
                        "✓ Hamiltonian\n"
                        "✓ Absorption & Fluorescence Spectra\n\n"
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
        self.pigment_labels = None

        self.analysis_dashboard.clear_all()
        self.data_table.clear_data()
        self.protein_viewer.clear_structure()
        self.hamiltonian_widget.clear()
        self.spectrum_widget.clear()
        self.hamiltonian_widget.set_run_enabled(False)
        self.spectrum_widget.set_run_enabled(False)

        self.status_message.emit("All data cleared")
