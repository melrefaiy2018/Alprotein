"""
Main Alprotein GUI Application
"""

import sys
import os
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                            QHBoxLayout, QSplitter, QTabWidget, QMenuBar, 
                            QStatusBar, QAction, QMessageBox, QProgressBar,
                            QLabel)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt5.QtGui import QIcon, QFont

from .main_window import MainWindow


class AlproteinGUI(QMainWindow):
    """
    Main Alprotein GUI Application Window
    
    This is the primary interface for the Alprotein package, providing
    an integrated environment for protein structure loading, site energy
    calculations, and result visualization.
    """
    
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Alprotein - Protein Site Energy Calculator")
        self.setGeometry(100, 100, 1400, 900)
        
        # Initialize the main window
        self.main_window = MainWindow()
        self.setCentralWidget(self.main_window)
        
        # Setup UI components
        self.setup_menubar()
        self.setup_statusbar()
        self.setup_window_properties()
        
        # Connect signals
        self.connect_signals()
        
        print("🧬 Alprotein GUI initialized successfully")
    
    def setup_menubar(self):
        """Setup the application menu bar"""
        menubar = self.menuBar()
        
        # File menu
        file_menu = menubar.addMenu('&File')
        
        # Load PDB action
        load_action = QAction('&Load PDB File...', self)
        load_action.setShortcut('Ctrl+O')
        load_action.setStatusTip('Load a PDB file for analysis')
        load_action.triggered.connect(self.main_window.load_pdb_file)
        file_menu.addAction(load_action)
        
        file_menu.addSeparator()
        
        # Export results action
        export_action = QAction('&Export Results...', self)
        export_action.setShortcut('Ctrl+E')
        export_action.setStatusTip('Export calculation results')
        export_action.triggered.connect(self.main_window.export_results)
        export_action.setEnabled(False)  # Enable when results are available
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
        site_energy_action.setShortcut('Ctrl+R')
        site_energy_action.setStatusTip('Calculate site energies for all pigments')
        site_energy_action.triggered.connect(self.main_window.calculate_site_energies)
        site_energy_action.setEnabled(False)  # Enable when structure is loaded
        calc_menu.addAction(site_energy_action)
        self.site_energy_action = site_energy_action
        
        # CDC analysis action
        cdc_action = QAction('&CDC Analysis', self)
        cdc_action.setShortcut('Ctrl+D')
        cdc_action.setStatusTip('Perform detailed contribution analysis')
        cdc_action.triggered.connect(self.main_window.run_cdc_analysis)
        cdc_action.setEnabled(False)  # Enable when results are available
        calc_menu.addAction(cdc_action)
        self.cdc_action = cdc_action
        
        # View menu
        view_menu = menubar.addMenu('&View')
        
        # Show Actions Window action
        show_actions_action = QAction('&Show Actions Window', self)
        show_actions_action.setStatusTip('Toggle visibility of the actions window')
        show_actions_action.triggered.connect(self.toggle_actions_window)
        view_menu.addAction(show_actions_action)
        
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
    
    def setup_statusbar(self):
        """Setup the status bar"""
        self.status_bar = self.statusBar()
        
        # Main status label
        self.status_label = QLabel("Ready")
        self.status_bar.addWidget(self.status_label)
        
        # Progress bar for calculations
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        self.progress_bar.setMaximumWidth(200)
        self.status_bar.addPermanentWidget(self.progress_bar)
        
        # Calculation status
        self.calc_status = QLabel("")
        self.status_bar.addPermanentWidget(self.calc_status)
    
    def setup_window_properties(self):
        """Setup window properties and styling"""
        # Set minimum size
        self.setMinimumSize(1000, 700)
        
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
                border-bottom: 1px solid #e1e1e1;
                padding: 6px 8px;
            }
            QMenuBar::item {
                color: #111111;
                padding: 6px 10px;
                margin: 2px;
            }
            QMenuBar::item:selected {
                background-color: #efefef;
                color: #111111;
                border-radius: 2px;
            }
            QMenu {
                background-color: #ffffff;
                color: #111111;
                border: 1px solid #e1e1e1;
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
            QLabel {
                color: #111111;
            }
        """)
    
    def connect_signals(self):
        """Connect signals between components"""
        # Connect main window signals to update menu states
        self.main_window.structure_loaded.connect(self.on_structure_loaded)
        self.main_window.calculation_started.connect(self.on_calculation_started)
        self.main_window.calculation_finished.connect(self.on_calculation_finished)
        self.main_window.status_message.connect(self.update_status)
        self.main_window.progress_update.connect(self.update_progress)
    
    def on_structure_loaded(self, success, message):
        """Handle structure loading completion"""
        if success:
            self.site_energy_action.setEnabled(True)
            self.update_status(f"Structure loaded: {message}")
        else:
            self.update_status(f"Failed to load structure: {message}")
            QMessageBox.warning(self, "Load Error", f"Failed to load structure:\n{message}")
    
    def on_calculation_started(self, calc_type):
        """Handle calculation start"""
        self.site_energy_action.setEnabled(False)
        self.cdc_action.setEnabled(False)
        self.export_action.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setValue(0)
        self.calc_status.setText(f"Calculating {calc_type}...")
        self.update_status(f"Running {calc_type} calculation...")
    
    def on_calculation_finished(self, calc_type, success, message):
        """Handle calculation completion"""
        self.progress_bar.setVisible(False)
        self.calc_status.setText("")
        
        if success:
            self.site_energy_action.setEnabled(True)
            self.cdc_action.setEnabled(True)
            self.export_action.setEnabled(True)
            self.update_status(f"{calc_type} calculation completed successfully")
            
            # Show completion message
            # QMessageBox.information(self, "Calculation Complete",  #TODO: t shoud be just a message
            #                       f"{calc_type} calculation completed successfully!\n\n{message}")
        else:
            self.site_energy_action.setEnabled(True)  # Re-enable for retry
            self.update_status(f"{calc_type} calculation failed")
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
        <h2>Alprotein</h2>
        <p><b>Protein Site Energy Calculator</b></p>
        <p>A comprehensive toolkit for calculating site energies and 
        analyzing pigment-protein interactions using the CDC method.</p>
        <p><b>Features:</b></p>
        <ul>
        <li>Site energy calculations using CDC method</li>
        <li>3D protein structure visualization</li>
        <li>Interactive contribution analysis</li>
        <li>Comprehensive result export</li>
        </ul>
        <p><i>Developed for research in photosynthetic protein complexes</i></p>
        """
        
        QMessageBox.about(self, "About Alprotein", about_text)
    
    def open_documentation(self):
        """Open documentation (placeholder)"""
        QMessageBox.information(
            self,
            "Documentation",
            "Documentation will be available at:\nhttps://github.com/melrefaiy2018/Alprotein-Alpha",
        )
    
    def toggle_actions_window(self):
        """Toggle the visibility of the action buttons window"""
        if self.main_window.action_buttons_window.isVisible():
            self.main_window.action_buttons_window.hide()
        else:
            self.main_window.action_buttons_window.show()
            self.main_window.action_buttons_window.raise_()
            self.main_window.action_buttons_window.activateWindow()
    
    def closeEvent(self, event):
        print("[DEBUG] closeEvent called: starting cleanup...")
        # Check if calculations are running
        if self.main_window.is_calculating():
            reply = QMessageBox.question(self, "Calculation in Progress", 
                                       "A calculation is currently running. Do you want to stop it and exit?",
                                       QMessageBox.Yes | QMessageBox.No)
            if reply == QMessageBox.No:
                event.ignore()
                return
            else:
                self.main_window.stop_calculation()
        
        # Save any necessary data
        self.main_window.cleanup()
        print("[DEBUG] Main window cleanup completed.")
        
        # Clear QtWebEngine cache and cookies
        try:
            from PyQt5.QtWebEngineWidgets import QWebEngineProfile
            web_view = getattr(self.main_window.protein_viewer.viewer_3d, 'web_view', None)
            if web_view:
                profile = web_view.page().profile()
                profile.clearHttpCache()
                profile.clearAllVisitedLinks()
                cookie_store = profile.cookieStore()
                if hasattr(cookie_store, 'deleteAllCookies'):
                    cookie_store.deleteAllCookies()
                else:
                    print("[INFO] QWebEngineCookieStore.deleteAllCookies() not available in this PyQt5 version.")
        except Exception as e:
            print(f"[INFO] Could not clear web engine cache: {e}")
        
        event.accept()
        from PyQt5.QtWidgets import QApplication
        QApplication.quit()
        import os
        os._exit(0)


def main():
    """Main entry point for the GUI application"""
    import os
    import sys
    
    # Suppress macOS NSOpenPanel warnings
    if sys.platform == 'darwin':
        os.environ['QT_MAC_WANTS_LAYER'] = '1'
    
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    
    # Configure matplotlib globally for black text
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
    app.setApplicationName("Alprotein")
    app.setApplicationVersion("1.0")
    app.setOrganizationName("Alprotein Team")
    
    # Create and show the main window
    gui = AlproteinGUI()
    gui.show()
    
    # Start the event loop
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
