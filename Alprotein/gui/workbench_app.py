"""
Scientific Workbench Application

Main application for the Alprotein Scientific Workbench with:
- Professional menu bar
- Status bar with progress tracking
- Keyboard shortcuts
- Help system
"""

import sys
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QAction, QMessageBox,
    QProgressBar, QLabel, QFileDialog
)
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QKeySequence, QIcon, QFont
from pathlib import Path

from Alprotein.gui.workbench_window import ScientificWorkbenchWindow
from Alprotein.gui.styles import GLOBAL_STYLESHEET, CHART_STYLE
import matplotlib.pyplot as plt

class ScientificWorkbenchApp(QMainWindow):
    """
    Main application window for Scientific Workbench
    """

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Alprotein Scientific Workbench")
        self.setGeometry(50, 50, 1600, 1000)

        # Central widget
        self.workbench = ScientificWorkbenchWindow()
        self.setCentralWidget(self.workbench)

        # Setup UI components
        self.setup_menu_bar()
        self.setup_status_bar()
        self.setup_shortcuts()
        self.connect_signals()

        # Apply styling
        self.apply_styling()

    def setup_menu_bar(self):
        """Create professional menu bar"""
        menubar = self.menuBar()
        menubar.setStyleSheet("""
            QMenuBar {
                background-color: #ffffff;
                color: #111111;
                padding: 6px 8px;
                font-size: 11px;
                border-bottom: 1px solid #e1e1e1;
            }
            QMenuBar::item {
                background-color: transparent;
                padding: 6px 10px;
            }
            QMenuBar::item:selected {
                background-color: #efefef;
            }
            QMenu {
                background-color: #ffffff;
                color: #111111;
                border: 1px solid #e1e1e1;
            }
            QMenu::item {
                padding: 6px 25px;
            }
            QMenu::item:selected {
                background-color: #111111;
                color: #ffffff;
            }
        """)

        # File Menu
        file_menu = menubar.addMenu("File")

        self.open_action = QAction("Open PDB...", self)
        self.open_action.setShortcut(QKeySequence.Open)
        self.open_action.setStatusTip("Open a PDB file")
        self.open_action.triggered.connect(self.workbench.on_open_project)
        file_menu.addAction(self.open_action)

        file_menu.addSeparator()

        self.export_action = QAction("Export Data...", self)
        self.export_action.setShortcut("Ctrl+E")
        self.export_action.setStatusTip("Export calculation results")
        self.export_action.triggered.connect(self.on_export_data)
        file_menu.addAction(self.export_action)

        self.export_image_action = QAction("Export Plots...", self)
        self.export_image_action.setShortcut("Ctrl+Shift+E")
        self.export_image_action.setStatusTip("Export plots as images")
        self.export_image_action.triggered.connect(self.on_export_plots)
        file_menu.addAction(self.export_image_action)

        file_menu.addSeparator()

        self.close_action = QAction("Close Project", self)
        self.close_action.setShortcut("Ctrl+W")
        self.close_action.setStatusTip("Close current project")
        self.close_action.triggered.connect(self.on_close_project)
        file_menu.addAction(self.close_action)

        file_menu.addSeparator()

        self.exit_action = QAction("Exit", self)
        self.exit_action.setShortcut(QKeySequence.Quit)
        self.exit_action.setStatusTip("Exit application")
        self.exit_action.triggered.connect(self.close)
        file_menu.addAction(self.exit_action)

        # Edit Menu
        edit_menu = menubar.addMenu("Edit")

        self.settings_action = QAction("Settings...", self)
        self.settings_action.setShortcut("Ctrl+,")
        self.settings_action.setStatusTip("Open settings")
        self.settings_action.triggered.connect(self.on_settings)
        edit_menu.addAction(self.settings_action)

        # Structure Menu
        structure_menu = menubar.addMenu("Structure")

        self.properties_action = QAction("Properties...", self)
        self.properties_action.setShortcut("Ctrl+I")
        self.properties_action.setStatusTip("View structure properties")
        self.properties_action.triggered.connect(self.workbench.on_show_properties)
        structure_menu.addAction(self.properties_action)

        structure_menu.addSeparator()

        self.validate_action = QAction("Validate Structure", self)
        self.validate_action.setShortcut("Ctrl+Shift+V")
        self.validate_action.setStatusTip("Validate current structure")
        self.validate_action.triggered.connect(self.on_validate_structure)
        structure_menu.addAction(self.validate_action)

        # Calculate Menu
        calc_menu = menubar.addMenu("Calculate")

        self.calc_site_energy_action = QAction("Site Energies", self)
        self.calc_site_energy_action.setShortcut("Ctrl+1")
        self.calc_site_energy_action.setStatusTip("Calculate site energies")
        self.calc_site_energy_action.triggered.connect(
            lambda: self.workbench.on_run_calculation("site_energies")
        )
        calc_menu.addAction(self.calc_site_energy_action)

        self.calc_hamiltonian_action = QAction("Hamiltonian", self)
        self.calc_hamiltonian_action.setShortcut("Ctrl+2")
        self.calc_hamiltonian_action.setStatusTip("Construct Hamiltonian matrix")
        self.calc_hamiltonian_action.triggered.connect(
            lambda: self.workbench.on_run_calculation("hamiltonian")
        )
        calc_menu.addAction(self.calc_hamiltonian_action)

        self.calc_spectrum_action = QAction("Absorption Spectrum", self)
        self.calc_spectrum_action.setShortcut("Ctrl+3")
        self.calc_spectrum_action.setStatusTip("Calculate absorption spectrum")
        self.calc_spectrum_action.triggered.connect(
            lambda: self.workbench.on_run_calculation("spectrum")
        )
        calc_menu.addAction(self.calc_spectrum_action)

        calc_menu.addSeparator()

        self.calc_all_action = QAction("Run All Calculations", self)
        self.calc_all_action.setShortcut("Ctrl+Shift+R")
        self.calc_all_action.setStatusTip("Run complete workflow")
        self.calc_all_action.triggered.connect(self.on_run_all_calculations)
        calc_menu.addAction(self.calc_all_action)

        # Analyze Menu
        analyze_menu = menubar.addMenu("Analyze")

        self.view_site_energies_action = QAction("Site Energy Analysis", self)
        self.view_site_energies_action.setShortcut("Ctrl+Shift+1")
        self.view_site_energies_action.setStatusTip("View detailed site energy analysis")
        self.view_site_energies_action.triggered.connect(
            lambda: self.workbench.on_view_details("site_energies")
        )
        analyze_menu.addAction(self.view_site_energies_action)

        self.view_hamiltonian_action = QAction("Hamiltonian Analysis", self)
        self.view_hamiltonian_action.setShortcut("Ctrl+Shift+2")
        self.view_hamiltonian_action.setStatusTip("View detailed Hamiltonian analysis")
        self.view_hamiltonian_action.triggered.connect(
            lambda: self.workbench.on_view_details("hamiltonian")
        )
        analyze_menu.addAction(self.view_hamiltonian_action)

        self.view_spectrum_action = QAction("Spectrum Analysis", self)
        self.view_spectrum_action.setShortcut("Ctrl+Shift+3")
        self.view_spectrum_action.setStatusTip("View detailed spectrum analysis")
        self.view_spectrum_action.triggered.connect(
            lambda: self.workbench.on_view_details("spectrum")
        )
        analyze_menu.addAction(self.view_spectrum_action)

        # Visualize Menu
        viz_menu = menubar.addMenu("Visualize")

        self.reset_view_action = QAction("Reset 3D View", self)
        self.reset_view_action.setShortcut("Ctrl+R")
        self.reset_view_action.setStatusTip("Reset 3D viewer to default view")
        self.reset_view_action.triggered.connect(self.on_reset_view)
        viz_menu.addAction(self.reset_view_action)

        viz_menu.addSeparator()

        # Tab switching actions
        self.goto_structure_action = QAction("Go to 3D Structure", self)
        self.goto_structure_action.setShortcut("Ctrl+Shift+1")
        self.goto_structure_action.setStatusTip("Switch to 3D Structure tab")
        self.goto_structure_action.triggered.connect(lambda: self.workbench.workspace_tabs.setCurrentIndex(0))
        viz_menu.addAction(self.goto_structure_action)

        self.goto_dashboard_action = QAction("Go to Analysis Dashboard", self)
        self.goto_dashboard_action.setShortcut("Ctrl+Shift+2")
        self.goto_dashboard_action.setStatusTip("Switch to Analysis Dashboard tab")
        self.goto_dashboard_action.triggered.connect(lambda: self.workbench.workspace_tabs.setCurrentIndex(1))
        viz_menu.addAction(self.goto_dashboard_action)

        self.goto_table_action = QAction("Go to Data Table", self)
        self.goto_table_action.setShortcut("Ctrl+Shift+3")
        self.goto_table_action.setStatusTip("Switch to Data Table tab")
        self.goto_table_action.triggered.connect(lambda: self.workbench.workspace_tabs.setCurrentIndex(2))
        viz_menu.addAction(self.goto_table_action)

        # Help Menu
        help_menu = menubar.addMenu("Help")

        self.quick_start_action = QAction("Quick Start Guide", self)
        self.quick_start_action.setShortcut(QKeySequence.HelpContents)
        self.quick_start_action.setStatusTip("View quick start guide")
        self.quick_start_action.triggered.connect(self.on_quick_start)
        help_menu.addAction(self.quick_start_action)

        help_menu.addSeparator()

        self.about_action = QAction("About Alprotein", self)
        self.about_action.setStatusTip("About this application")
        self.about_action.triggered.connect(self.on_about)
        help_menu.addAction(self.about_action)

    def setup_status_bar(self):
        """Create status bar with progress tracking"""
        statusbar = self.statusBar()
        statusbar.setStyleSheet("""
            QStatusBar {
                background-color: #ffffff;
                color: #111111;
                border-top: 1px solid #e1e1e1;
                font-size: 10px;
            }
        """)

        # Status label
        self.status_label = QLabel("Ready")
        self.status_label.setStyleSheet("padding: 4px 8px;")
        statusbar.addWidget(self.status_label, 1)

        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setMaximumWidth(200)
        self.progress_bar.setMaximumHeight(16)
        self.progress_bar.setVisible(False)
        self.progress_bar.setStyleSheet("""
            QProgressBar {
                border: 1px solid #e1e1e1;
                border-radius: 2px;
                text-align: center;
                background-color: #f3f3f3;
            }
            QProgressBar::chunk {
                background-color: #111111;
            }
        """)
        statusbar.addPermanentWidget(self.progress_bar)

        # Memory/info label
        self.info_label = QLabel("")
        self.info_label.setStyleSheet("padding: 4px 8px; color: #6b6b6b;")
        statusbar.addPermanentWidget(self.info_label)

    def setup_shortcuts(self):
        """Setup additional keyboard shortcuts"""
        # Already setup in menu actions
        pass

    def connect_signals(self):
        """Connect workbench signals to app"""
        self.workbench.status_message.connect(self.update_status)
        self.workbench.progress_update.connect(self.update_progress)

    def apply_styling(self):
        """Apply global application styling"""
        # Apply chart styling
        for key, value in CHART_STYLE.items():
            plt.rcParams[key] = value

    def update_status(self, message: str):
        """Update status bar message"""
        self.status_label.setText(message)

        # Auto-clear after 5 seconds for non-error messages
        if not message.startswith("Error"):
            QTimer.singleShot(5000, lambda: self.status_label.setText("Ready"))

    def update_progress(self, message: str, percentage: int):
        """Update progress bar"""
        if percentage > 0 and percentage < 100:
            self.progress_bar.setVisible(True)
            self.progress_bar.setValue(percentage)
            self.status_label.setText(message)
        else:
            self.progress_bar.setVisible(False)

    def on_export_data(self):
        """Export calculation data"""
        self.workbench.on_export_table()

    def on_export_plots(self):
        """Export plots as images"""
        QMessageBox.information(
            self,
            "Export Plots",
            "Click the 'Export' button on each plot in the Analysis Dashboard to export individual plots."
        )

    def on_close_project(self):
        """Close current project"""
        reply = QMessageBox.question(
            self,
            "Close Project",
            "Are you sure you want to close the current project?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )

        if reply == QMessageBox.Yes:
            self.workbench.clear_all_data()
            self.update_status("Project closed")

    def on_settings(self):
        """Open settings dialog"""
        QMessageBox.information(
            self,
            "Settings",
            "Settings are available in the Tools Panel on the right.\n\n"
            "Adjust parameters in the SETTINGS section:\n"
            "- Dielectric constants\n"
            "- Oscillator strength\n"
            "- Temperature\n"
            "- Coupling method"
        )

    def on_validate_structure(self):
        """Validate current structure"""
        if not self.workbench.pigment_system:
            QMessageBox.warning(
                self,
                "No Structure",
                "Please load a PDB file first"
            )
            return

        # Structure is validated on load, so just show info
        self.workbench.on_show_properties()

    def on_run_all_calculations(self):
        """Run complete calculation workflow"""
        if not self.workbench.pigment_system:
            QMessageBox.warning(
                self,
                "No Structure",
                "Please load a PDB file first"
            )
            return

        reply = QMessageBox.question(
            self,
            "Run All Calculations",
            "This will run the complete workflow:\n"
            "1. Site Energies\n"
            "2. Hamiltonian Construction\n"
            "3. Spectrum Calculation\n\n"
            "This may take several minutes. Continue?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )

        if reply == QMessageBox.Yes:
            # Start with site energies
            # The workflow will continue automatically via signals
            self.workbench.on_run_calculation("site_energies")

    def on_reset_view(self):
        """Reset 3D view"""
        self.workbench.protein_viewer.reset_view()
        self.update_status("3D view reset")

    def on_quick_start(self):
        """Show quick start guide"""
        guide = """
        <h2>Alprotein Scientific Workbench - Quick Start</h2>

        <h3>1. Load Structure</h3>
        <p>• File → Open PDB... (or Ctrl+O)<br>
        • Select your PDB file<br>
        • View structure in 3D viewer (first tab)</p>

        <h3>2. Configure Parameters</h3>
        <p>• Adjust settings in Tools Panel (right sidebar)<br>
        • Set dielectric constants<br>
        • Choose coupling method (TrEsp or Dipole)</p>

        <h3>3. Run Calculations</h3>
        <p>• <b>Site Energies</b>: Calculate → Site Energies (Ctrl+1)<br>
        • <b>Hamiltonian</b>: Calculate → Hamiltonian (Ctrl+2)<br>
        • <b>Spectrum</b>: Calculate → Absorption Spectrum (Ctrl+3)</p>

        <h3>4. View Results (Tabbed Interface)</h3>
        <p>• <b>Tab 1 - 3D Structure</b>: Protein visualization (Ctrl+Shift+1)<br>
        • <b>Tab 2 - Analysis Dashboard</b>: Three mini-plots (Ctrl+Shift+2)<br>
        • <b>Tab 3 - Data Table</b>: Detailed data (Ctrl+Shift+3)<br>
        • Click "Details" or "Export" in dashboard for full views</p>

        <h3>Keyboard Shortcuts</h3>
        <p><b>File Operations:</b><br>
        • <b>Ctrl+O</b>: Open PDB<br>
        • <b>Ctrl+E</b>: Export data<br>
        <b>Calculations:</b><br>
        • <b>Ctrl+1/2/3</b>: Run site energies/Hamiltonian/spectrum<br>
        <b>View Tabs:</b><br>
        • <b>Ctrl+Shift+1/2/3</b>: Switch to Structure/Dashboard/Table<br>
        <b>Other:</b><br>
        • <b>Ctrl+R</b>: Reset 3D view<br>
        • <b>F1</b>: This help</p>

        <h3>Tips</h3>
        <p>• Switch between tabs to view different aspects of your analysis<br>
        • Watch progress in status bar during calculations<br>
        • Select pigments in table to highlight in 3D viewer<br>
        • Use Ctrl+Shift+1/2/3 for quick tab switching</p>
        """

        msg = QMessageBox(self)
        msg.setWindowTitle("Quick Start Guide")
        msg.setTextFormat(Qt.RichText)
        msg.setText(guide)
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()

    def on_about(self):
        """Show about dialog"""
        about_text = """
        <h2>Alprotein Scientific Workbench</h2>
        <p><b>Version 2.0</b></p>

        <p>Professional software for pigment-protein complex spectroscopy analysis.</p>

        <h3>Features</h3>
        <ul>
        <li>Site energy calculation using CDC method</li>
        <li>Hamiltonian construction (TrEsp/Dipole coupling)</li>
        <li>Linear absorption spectrum calculation</li>
        <li>3D molecular visualization</li>
        <li>Professional data analysis and export</li>
        </ul>

        <h3>Components</h3>
        <p>• 3D Structure Viewer<br>
        • Analysis Dashboard<br>
        • Data Table<br>
        • Tools Panel</p>

        <p><small>Built with PyQt5 and Matplotlib</small></p>
        """

        msg = QMessageBox(self)
        msg.setWindowTitle("About Alprotein")
        msg.setTextFormat(Qt.RichText)
        msg.setText(about_text)
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()


def main():
    """Main application entry point"""
    app = QApplication(sys.argv)

    # Set application metadata
    app.setApplicationName("Alprotein Scientific Workbench")
    app.setOrganizationName("Alprotein")
    app.setApplicationVersion("2.0")

    # Set default font
    font = QFont("Segoe UI", 10)
    app.setFont(font)

    # Apply global stylesheet
    app.setStyleSheet(GLOBAL_STYLESHEET)

    # Create and show main window
    window = ScientificWorkbenchApp()
    window.show()

    # Start event loop
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
