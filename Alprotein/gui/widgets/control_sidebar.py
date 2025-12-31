"""
Consolidated and Tabbed Control Sidebar Widget
"""

from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QGroupBox, QSizePolicy, QPushButton, QHBoxLayout, QLabel, QScrollArea)
from PyQt5.QtCore import pyqtSignal, Qt
from PyQt5.QtGui import QFont

from .file_loader import FileLoaderWidget
from .calculation_panel import CalculationPanel
from .control_panel_3d import ControlPanel3D
from .actions_panel import ActionsPanel

class CollapsibleBox(QWidget):
    """A collapsible group box widget."""
    def __init__(self, title="", parent=None):
        super(CollapsibleBox, self).__init__(parent)
        self.toggle_button = QPushButton(title)
        self.toggle_button.setStyleSheet("""
            QPushButton {
                text-align: left;
                font-weight: bold;
                border: none;
                padding: 5px;
            }
            QPushButton:checked {
                background-color: #efefef;
            }
        """)
        self.toggle_button.setCheckable(True)
        self.toggle_button.setChecked(True)

        self.content_area = QWidget()
        self.content_layout = QVBoxLayout(self.content_area)
        
        main_layout = QVBoxLayout(self)
        main_layout.setSpacing(0)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.addWidget(self.toggle_button)
        main_layout.addWidget(self.content_area)

        self.toggle_button.toggled.connect(self.toggle)

    def setContentLayout(self, layout):
        # Clear existing layout
        while self.content_layout.count():
            item = self.content_layout.takeAt(0)
            widget = item.widget()
            if widget:
                widget.setParent(None)
        
        # Add new layout
        self.content_layout.addLayout(layout)

    def toggle(self, checked):
        self.content_area.setVisible(checked)

class ControlSidebar(QWidget):
    """
    A redesigned sidebar with collapsible sections for a modern, clean look.
    """
    # Signals to be exposed to MainWindow
    file_selected = pyqtSignal(str)
    calculate_site_energies = pyqtSignal()
    run_cdc_analysis = pyqtSignal()
    settings_changed = pyqtSignal(dict)
    view_changed = pyqtSignal(dict)

    def __init__(self):
        super().__init__()
        self.setFixedWidth(450) # control the default width of the side bar on left 
        self.setup_ui()
        self.connect_signals()

    def setup_ui(self):
        """Setup the redesigned control sidebar UI."""
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll.setStyleSheet("QScrollBar:vertical { width: 8px; margin: 2px; }")

        content_widget = QWidget()
        content_layout = QVBoxLayout(content_widget)
        content_layout.setSpacing(15)
        content_layout.setContentsMargins(16, 16, 16, 16)
        content_layout.setAlignment(Qt.AlignTop)

        # --- Load PDB Section ---
        self.load_pdb_box = CollapsibleBox("Load PDB")
        load_pdb_layout = QVBoxLayout()
        self.file_loader = FileLoaderWidget()
        load_pdb_layout.addWidget(self.file_loader)
        self.load_pdb_box.setContentLayout(load_pdb_layout)
        content_layout.addWidget(self.load_pdb_box)

        # --- System Status Section ---
        self.system_status_box = CollapsibleBox("System Status")
        system_status_layout = QVBoxLayout()
        self.status_label = QLabel("No PDB loaded.")
        self.status_label.setFont(QFont("Helvetica", 10))
        system_status_layout.addWidget(self.status_label)
        self.system_status_box.setContentLayout(system_status_layout)
        content_layout.addWidget(self.system_status_box)

        # --- Calculation Settings Section ---
        self.calc_settings_box = CollapsibleBox("Calculation Settings")
        calc_settings_layout = QVBoxLayout()
        self.calculation_panel = CalculationPanel()
        calc_settings_layout.addWidget(self.calculation_panel)
        self.calc_settings_box.setContentLayout(calc_settings_layout)
        content_layout.addWidget(self.calc_settings_box)

        # --- Actions Section ---
        self.actions_box = CollapsibleBox("Actions")
        actions_layout = QVBoxLayout()
        self.actions_panel = ActionsPanel()
        actions_layout.addWidget(self.actions_panel)
        self.actions_box.setContentLayout(actions_layout)
        content_layout.addWidget(self.actions_box)

        # --- Advanced Options Section ---
        self.advanced_options_box = CollapsibleBox("Advanced Options")
        advanced_options_layout = QVBoxLayout()
        self.control_panel_3d = ControlPanel3D()
        advanced_options_layout.addWidget(self.control_panel_3d)
        self.advanced_options_box.setContentLayout(advanced_options_layout)
        content_layout.addWidget(self.advanced_options_box)

        content_layout.addStretch(1)
        scroll.setWidget(content_widget)
        main_layout.addWidget(scroll)

    def connect_signals(self):
        """Connect signals from sub-widgets to expose them."""
        self.file_loader.file_selected.connect(self.file_selected.emit)
        self.actions_panel.calculate_site_energies.connect(self.calculate_site_energies.emit)
        self.actions_panel.run_cdc_analysis.connect(self.run_cdc_analysis.emit)
        self.calculation_panel.settings_changed.connect(self.settings_changed.emit)
        self.control_panel_3d.view_changed.connect(self.view_changed.emit)

    def set_loaded_file(self, file_path: str, pigment_count: int):
        """Update UI elements with loaded file information."""
        self.file_loader.set_loaded_file(file_path, pigment_count)
        self.actions_panel.set_system_loaded(True)
        self.status_label.setText(f"{pigment_count} pigments loaded.")

    def set_system_info(self, validation: dict):
        """Update system status with validation info."""
        warnings = validation.get('warnings', [])
        if warnings:
            self.status_label.setText(f"System has warnings: {', '.join(warnings)}")
        else:
            self.status_label.setText("System is valid.")

    def update_pigment_list(self, pigments: dict):
        """Pass through to 3D control panel."""
        self.control_panel_3d.update_pigment_list(pigments)

    def set_site_energies_calculated(self, calculated: bool):
        """Pass through to actions panel."""
        self.actions_panel.set_site_energies_calculated(calculated)

    def get_calculation_settings(self) -> dict:
        """Get current calculation settings from the panel."""
        return self.calculation_panel.get_current_settings()
    
    def set_calculation_running(self, running: bool, calc_type: str = ""):
        """Pass through to actions panel."""
        self.actions_panel.set_calculation_running(running, calc_type)
    
    def update_progress(self, percentage: int, message: str):
        """Pass through to actions panel."""
        self.actions_panel.update_progress(percentage, message)
