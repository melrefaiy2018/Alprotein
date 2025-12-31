"""
Alprotein GUI Module

This module provides a graphical user interface for protein site energy calculations
and visualization using the Alprotein package.
"""

from .app import AlproteinGUI
from .main_window import MainWindow

__all__ = ['AlproteinGUI', 'MainWindow', 'launch_gui']


def launch_gui():
    """
    Launch the Alprotein GUI application.
    
    This is the main entry point for the GUI interface.
    """
    import sys
    from PyQt5.QtWidgets import QApplication
    
    # Create QApplication if it doesn't exist
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    
    # Create and show the main window
    gui = AlproteinGUI()
    gui.show()
    
    # Start the event loop and exit properly
    sys.exit(app.exec_())
