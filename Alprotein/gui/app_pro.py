"""
Professional Alprotein GUI Application

This is the main application window with the enhanced professional layout.
"""

import sys
import os
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout,
    QMenuBar, QStatusBar, QAction, QMessageBox, QProgressBar,
    QLabel
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont, QIcon

from .main_window_pro import MainWindowPro


class AlproteinGUIPro(QMainWindow):
    """
    Professional Alprotein GUI Application Window
    """

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Alprotein - Pigment-Protein Complex Spectroscopy")
        self.setGeometry(50, 50, 1600, 1000)

        # Initialize the main window
        self.main_window = MainWindowPro()
        self.setCentralWidget(self.main_window)

        # Setup UI components
        self.setup_menubar()
        self.setup_statusbar()
        self.setup_window_properties()

        # Connect signals
        self.connect_signals()

        print("🎨 Alprotein Professional GUI initialized successfully")

    def setup_menubar(self):
        """Setup the application menu bar"""
        menubar = self.menuBar()

        # File menu
        file_menu = menubar.addMenu('&File')

        # Load PDB action
        load_action = QAction('&Load PDB File...', self)
        load_action.setShortcut('Ctrl+O')
        load_action.setStatusTip('Load a PDB file for analysis')
        load_action.triggered.connect(self.load_pdb_file)
        file_menu.addAction(load_action)

        file_menu.addSeparator()

        # Export results action
        export_action = QAction('&Export Results...', self)
        export_action.setShortcut('Ctrl+E')
        export_action.setStatusTip('Export calculation results')
        export_action.triggered.connect(self.export_results)
        file_menu.addAction(export_action)
        self.export_action = export_action

        file_menu.addSeparator()

        # Exit action
        exit_action = QAction('E&xit', self)
        exit_action.setShortcut('Ctrl+Q')
        exit_action.setStatusTip('Exit the application')
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        # Calculate menu
        calc_menu = menubar.addMenu('&Calculate')

        # Site energies action
        site_energy_action = QAction('&Site Energies', self)
        site_energy_action.setShortcut('Ctrl+1')
        site_energy_action.triggered.connect(self.main_window.calculate_site_energies)
        calc_menu.addAction(site_energy_action)

        # Hamiltonian action
        hamiltonian_action = QAction('&Hamiltonian', self)
        hamiltonian_action.setShortcut('Ctrl+2')
        hamiltonian_action.triggered.connect(self.main_window.construct_hamiltonian)
        calc_menu.addAction(hamiltonian_action)

        # Diagonalize action
        diag_action = QAction('&Diagonalize', self)
        diag_action.setShortcut('Ctrl+3')
        diag_action.triggered.connect(self.main_window.diagonalize_hamiltonian)
        calc_menu.addAction(diag_action)

        # Spectrum action
        spectrum_action = QAction('&Absorption Spectrum', self)
        spectrum_action.setShortcut('Ctrl+4')
        spectrum_action.triggered.connect(self.main_window.calculate_spectrum)
        calc_menu.addAction(spectrum_action)

        # View menu
        view_menu = menubar.addMenu('&View')

        # Zoom actions
        zoom_in_action = QAction('Zoom &In', self)
        zoom_in_action.setShortcut('Ctrl++')
        view_menu.addAction(zoom_in_action)

        zoom_out_action = QAction('Zoom &Out', self)
        zoom_out_action.setShortcut('Ctrl+-')
        view_menu.addAction(zoom_out_action)

        reset_view_action = QAction('&Reset View', self)
        reset_view_action.setShortcut('Ctrl+0')
        view_menu.addAction(reset_view_action)

        # Help menu
        help_menu = menubar.addMenu('&Help')

        # About action
        about_action = QAction('&About Alprotein', self)
        about_action.setStatusTip('About this application')
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)

        # Documentation action
        docs_action = QAction('&Documentation', self)
        docs_action.setStatusTip('Open documentation')
        docs_action.triggered.connect(self.open_documentation)
        help_menu.addAction(docs_action)

        # Quick Start Guide
        quickstart_action = QAction('&Quick Start Guide', self)
        quickstart_action.triggered.connect(self.show_quickstart)
        help_menu.addAction(quickstart_action)

    def setup_statusbar(self):
        """Setup the status bar"""
        self.status_bar = self.statusBar()

        # Main status label
        self.status_label = QLabel("Ready")
        self.status_bar.addWidget(self.status_label)

        # Progress bar for calculations
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        self.progress_bar.setMaximumWidth(250)
        self.progress_bar.setStyleSheet("""
            QProgressBar {
                border: 1px solid #e1e1e1;
                border-radius: 2px;
                text-align: center;
                background-color: #f3f3f3;
                color: #111111;
            }
            QProgressBar::chunk {
                background-color: #111111;
            }
        """)
        self.status_bar.addPermanentWidget(self.progress_bar)

        # Calculation status
        self.calc_status = QLabel("")
        self.calc_status.setStyleSheet("color: #1f7a44; font-weight: bold;")
        self.status_bar.addPermanentWidget(self.calc_status)

    def setup_window_properties(self):
        """Setup window properties and styling"""
        # Set minimum size
        self.setMinimumSize(1200, 800)

        # Apply styling
        self.setStyleSheet("""
            QMainWindow {
                background-color: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                                                  stop:0 #fafafa, stop:1 #f2f2f2);
                color: #111111;
            }
            QMenuBar {
                background-color: #ffffff;
                color: #111111;
                padding: 6px 8px;
                border-bottom: 1px solid #e1e1e1;
            }
            QMenuBar::item {
                color: #111111;
                padding: 6px 10px;
                margin: 0px 2px;
                background-color: transparent;
            }
            QMenuBar::item:selected {
                background-color: #efefef;
                border-radius: 2px;
            }
            QMenuBar::item:pressed {
                background-color: #111111;
                color: #ffffff;
            }
            QMenu {
                background-color: #ffffff;
                color: #111111;
                border: 1px solid #e1e1e1;
            }
            QMenu::item {
                padding: 6px 25px;
                color: #111111;
            }
            QMenu::item:selected {
                background-color: #111111;
                color: #ffffff;
            }
            QStatusBar {
                background-color: #ffffff;
                color: #111111;
                border-top: 1px solid #e1e1e1;
            }
            QStatusBar QLabel {
                color: #111111;
            }
        """)

    def connect_signals(self):
        """Connect signals between components"""
        self.main_window.structure_loaded.connect(self.on_structure_loaded)
        self.main_window.calculation_started.connect(self.on_calculation_started)
        self.main_window.calculation_finished.connect(self.on_calculation_finished)
        self.main_window.status_message.connect(self.update_status)
        self.main_window.progress_update.connect(self.update_progress)

    def load_pdb_file(self):
        """Open file dialog to load PDB file"""
        from PyQt5.QtWidgets import QFileDialog

        options = QFileDialog.DontUseNativeDialog if sys.platform == 'darwin' else QFileDialog.Options()

        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Load PDB File",
            "",
            "PDB Files (*.pdb);;All Files (*)",
            options=options
        )

        if file_path:
            self.main_window.load_pdb_file_path(file_path)

    def export_results(self):
        """Export calculation results"""
        QMessageBox.information(
            self,
            "Export Results",
            "Use the export buttons in each results tab to export specific data.\n\n"
            "Available exports:\n"
            "• Site Energies: CSV and plots\n"
            "• Hamiltonian: Matrix data\n"
            "• Eigenvalues: Eigenvectors CSV\n"
            "• Spectrum: Data and plots"
        )

    def on_structure_loaded(self, success, message):
        """Handle structure loading completion"""
        if success:
            self.update_status(f"✓ Structure loaded: {message}")
        else:
            self.update_status(f"✗ Failed to load structure")
            QMessageBox.warning(self, "Load Error", f"Failed to load structure:\n{message}")

    def on_calculation_started(self, calc_type):
        """Handle calculation start"""
        self.progress_bar.setVisible(True)
        self.progress_bar.setValue(0)
        self.calc_status.setText(f"Calculating {calc_type}...")
        self.update_status(f"▶ Running {calc_type} calculation...")

    def on_calculation_finished(self, calc_type, success, message):
        """Handle calculation completion"""
        self.progress_bar.setVisible(False)
        self.calc_status.setText("")

        if success:
            self.update_status(f"✓ {calc_type} calculation completed")
        else:
            self.update_status(f"✗ {calc_type} calculation failed")
            QMessageBox.critical(self, "Calculation Error",
                               f"{calc_type} calculation failed:\n{message}")

    def update_status(self, message):
        """Update the status bar message"""
        self.status_label.setText(message)
        print(f"Status: {message}")

    def update_progress(self, value, message=""):
        """Update the progress bar"""
        self.progress_bar.setValue(value)
        if message:
            self.calc_status.setText(message)

    def show_about(self):
        """Show about dialog"""
        about_text = """
        <div style="text-align: center;">
        <h2 style="color: #111111;">Alprotein</h2>
        <p style="font-size: 12px; color: #6b6b6b;"><b>Version 0.2.0</b></p>
        <hr>
        <p><b>Pigment-Protein Complex Spectroscopy Analysis</b></p>
        <p>A comprehensive platform for calculating excitonic properties and
        optical spectra of pigment-protein complexes.</p>

        <h3 style="color: #111111;">Features</h3>
        <ul style="text-align: left;">
        <li>Site energy calculations (CDC method)</li>
        <li>Hamiltonian construction & diagonalization</li>
        <li>Linear absorption spectra calculation</li>
        <li>3D protein structure visualization</li>
        <li>Interactive contribution analysis</li>
        <li>Professional results export</li>
        </ul>

        <hr>
        <p style="font-size: 10px; color: #8a8a8a;">
        Developed for research in photosynthetic protein complexes<br>
        © 2025 Alprotein Team
        </p>
        </div>
        """

        QMessageBox.about(self, "About Alprotein", about_text)

    def show_quickstart(self):
        """Show quick start guide"""
        quickstart_text = """
        <h3>Quick Start Guide</h3>

        <h4>1. Load Structure</h4>
        <p>• Click "Browse" in the Workflow Panel<br>
        • Select a PDB file containing pigment residues (CLA, CHL, etc.)</p>

        <h4>2. Calculate Site Energies</h4>
        <p>• Click "Calculate" in Step 3 of the Workflow Panel<br>
        • View results in the "Site Energy" tab</p>

        <h4>3. Construct Hamiltonian</h4>
        <p>• Click "Calculate" in Step 4<br>
        • View the matrix heatmap in "Hamiltonian" tab</p>

        <h4>4. Diagonalize & Get Spectrum</h4>
        <p>• Click "Diagonalize" to compute eigenvalues<br>
        • Click "Calculate" in Step 5 for absorption spectrum</p>

        <h4>5. Export Results</h4>
        <p>• Use export buttons in each results tab<br>
        • Available formats: CSV, PNG, PDF</p>

        <hr>
        <p><b>Keyboard Shortcuts:</b><br>
        Ctrl+O: Load PDB file<br>
        Ctrl+1-4: Run calculations<br>
        Ctrl+E: Export results</p>
        """

        msg = QMessageBox(self)
        msg.setWindowTitle("Quick Start Guide")
        msg.setTextFormat(Qt.RichText)
        msg.setText(quickstart_text)
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()

    def open_documentation(self):
        """Open documentation"""
        QMessageBox.information(
            self,
            "Documentation",
            "Full documentation available at:\n\n"
            "https://github.com/melrefaiy2018/Alprotein-Alpha\n\n"
            "For quick help, see Help → Quick Start Guide"
        )

    def closeEvent(self, event):
        """Handle window close event"""
        print("[DEBUG] closeEvent called: starting cleanup...")

        if self.main_window.is_calculating():
            reply = QMessageBox.question(
                self,
                "Calculation in Progress",
                "A calculation is currently running. Do you want to stop it and exit?",
                QMessageBox.Yes | QMessageBox.No
            )
            if reply == QMessageBox.No:
                event.ignore()
                return
            else:
                self.main_window.stop_calculation()

        # Cleanup
        self.main_window.cleanup()
        print("[DEBUG] Main window cleanup completed.")

        event.accept()
        QApplication.quit()


def main():
    """Main entry point for the professional GUI application"""
    import os
    import sys

    # Suppress macOS warnings
    if sys.platform == 'darwin':
        os.environ['QT_MAC_WANTS_LAYER'] = '1'

    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)

    # Configure matplotlib globally
    import matplotlib.pyplot as plt
    plt.style.use('default')
    plt.rcParams['text.color'] = 'black'
    plt.rcParams['axes.labelcolor'] = 'black'
    plt.rcParams['xtick.color'] = 'black'
    plt.rcParams['ytick.color'] = 'black'
    plt.rcParams['axes.edgecolor'] = 'black'
    plt.rcParams['figure.facecolor'] = 'white'
    plt.rcParams['axes.facecolor'] = 'white'

    # Set application properties
    app.setApplicationName("Alprotein Professional")
    app.setApplicationVersion("0.2.0")
    app.setOrganizationName("Alprotein Team")

    # Create and show the main window
    gui = AlproteinGUIPro()
    gui.show()

    # Start the event loop
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
