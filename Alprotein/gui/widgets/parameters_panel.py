"""
Parameters Panel Widget for Alprotein GUI

This widget provides a comprehensive interface for configuring calculation parameters:
- Method selection (TrEsp, Dipole-Dipole)
- Physical parameters (dielectric constants, oscillator strength, temperature)
- Pigment configuration table
"""

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QRadioButton, QDoubleSpinBox, QSpinBox,
    QTableWidget, QTableWidgetItem, QHeaderView, QComboBox,
    QFrame, QSizePolicy
)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont
from typing import Dict, Any


class ParametersPanel(QWidget):
    """
    Panel for configuring calculation parameters
    """

    # Signals
    parameters_changed = pyqtSignal(dict)  # parameter_dict
    pigment_added = pyqtSignal(dict)  # pigment_config
    pigment_removed = pyqtSignal(int)  # row_index

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setup_ui()

    def setup_ui(self):
        """Setup the UI"""
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(15, 15, 15, 15)
        main_layout.setSpacing(15)

        # Header
        header = QLabel("CALCULATION PARAMETERS")
        header.setAlignment(Qt.AlignCenter)
        header.setStyleSheet("""
            QLabel {
                background-color: #f3f3f3;
                color: #111111;
                padding: 10px;
                font-size: 13px;
                font-weight: bold;
                letter-spacing: 1px;
                border-radius: 2px;
                border: 1px solid #e1e1e1;
            }
        """)
        main_layout.addWidget(header)

        # Create horizontal layout for method and physical parameters
        top_layout = QHBoxLayout()

        # Method Selection Group
        method_group = self.create_method_group()
        top_layout.addWidget(method_group)

        # Physical Parameters Group
        physics_group = self.create_physics_group()
        top_layout.addWidget(physics_group)

        main_layout.addLayout(top_layout)

        # Pigment Configuration Group
        pigment_group = self.create_pigment_group()
        main_layout.addWidget(pigment_group, stretch=1)

        # Apply overall styling
        self.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                border: 1px solid #e1e1e1;
                border-radius: 2px;
                margin-top: 12px;
                padding-top: 10px;
                background-color: #ffffff;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px;
                color: #111111;
            }
            QLabel {
                color: #111111;
            }
            QDoubleSpinBox, QSpinBox, QComboBox {
                padding: 4px;
                border: 1px solid #e1e1e1;
                border-radius: 2px;
                background-color: #ffffff;
                color: #111111;
            }
            QDoubleSpinBox:focus, QSpinBox:focus, QComboBox:focus {
                border: 2px solid #111111;
            }
        """)

    def create_method_group(self):
        """Create method selection group"""
        group = QGroupBox("Method Selection")
        layout = QVBoxLayout(group)
        layout.setSpacing(8)

        # TrEsp method
        self.tresp_radio = QRadioButton("TrEsp (Transition Charge)")
        self.tresp_radio.setChecked(True)
        self.tresp_radio.setStyleSheet("color: black;")
        self.tresp_radio.toggled.connect(self.on_parameters_changed)
        layout.addWidget(self.tresp_radio)

        # Dipole-Dipole method
        self.dipole_radio = QRadioButton("Dipole-Dipole")
        self.dipole_radio.setStyleSheet("color: black;")
        self.dipole_radio.toggled.connect(self.on_parameters_changed)
        layout.addWidget(self.dipole_radio)

        layout.addStretch()

        return group

    def create_physics_group(self):
        """Create physical parameters group"""
        group = QGroupBox("Physical Parameters")
        layout = QVBoxLayout(group)
        layout.setSpacing(10)

        # Dielectric constant (CDC)
        cdc_layout = QHBoxLayout()
        cdc_layout.addWidget(QLabel("Dielectric (CDC):"))
        self.dielectric_cdc_spin = QDoubleSpinBox()
        self.dielectric_cdc_spin.setRange(1.0, 80.0)
        self.dielectric_cdc_spin.setValue(2.0)
        self.dielectric_cdc_spin.setDecimals(2)
        self.dielectric_cdc_spin.setSingleStep(0.1)
        self.dielectric_cdc_spin.valueChanged.connect(self.on_parameters_changed)
        cdc_layout.addWidget(self.dielectric_cdc_spin)
        layout.addLayout(cdc_layout)

        # Dielectric constant (TrEsp)
        tresp_layout = QHBoxLayout()
        tresp_layout.addWidget(QLabel("Dielectric (TrEsp):"))
        self.dielectric_tresp_spin = QDoubleSpinBox()
        self.dielectric_tresp_spin.setRange(1.0, 80.0)
        self.dielectric_tresp_spin.setValue(1.0)
        self.dielectric_tresp_spin.setDecimals(2)
        self.dielectric_tresp_spin.setSingleStep(0.1)
        self.dielectric_tresp_spin.valueChanged.connect(self.on_parameters_changed)
        tresp_layout.addWidget(self.dielectric_tresp_spin)
        layout.addLayout(tresp_layout)

        # Oscillator strength
        osc_layout = QHBoxLayout()
        osc_layout.addWidget(QLabel("Oscillator Strength:"))
        self.osc_strength_spin = QDoubleSpinBox()
        self.osc_strength_spin.setRange(0.1, 2.0)
        self.osc_strength_spin.setValue(0.72)
        self.osc_strength_spin.setDecimals(3)
        self.osc_strength_spin.setSingleStep(0.01)
        self.osc_strength_spin.valueChanged.connect(self.on_parameters_changed)
        osc_layout.addWidget(self.osc_strength_spin)
        layout.addLayout(osc_layout)

        # Temperature
        temp_layout = QHBoxLayout()
        temp_layout.addWidget(QLabel("Temperature (K):"))
        self.temperature_spin = QSpinBox()
        self.temperature_spin.setRange(0, 500)
        self.temperature_spin.setValue(300)
        self.temperature_spin.setSingleStep(10)
        self.temperature_spin.valueChanged.connect(self.on_parameters_changed)
        temp_layout.addWidget(self.temperature_spin)
        layout.addLayout(temp_layout)

        layout.addStretch()

        return group

    def create_pigment_group(self):
        """Create pigment configuration group"""
        group = QGroupBox("Pigment Configuration")
        layout = QVBoxLayout(group)
        layout.setSpacing(10)

        # Pigment table
        self.pigment_table = QTableWidget()
        self.pigment_table.setColumnCount(5)
        self.pigment_table.setHorizontalHeaderLabels([
            "Residue", "Pigment", "TrEsp Dict", "CDC Dict", "E_vac (cm⁻¹)"
        ])

        # Set column widths
        header = self.pigment_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.Stretch)
        header.setSectionResizeMode(2, QHeaderView.Stretch)
        header.setSectionResizeMode(3, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(4, QHeaderView.ResizeToContents)

        # Set alternating row colors
        self.pigment_table.setAlternatingRowColors(True)
        self.pigment_table.setStyleSheet("""
            QTableWidget {
                background-color: #ffffff;
                alternate-background-color: #f7f7f7;
                gridline-color: #e6e6e6;
                color: #111111;
            }
            QHeaderView::section {
                background-color: #f3f3f3;
                color: #111111;
                padding: 6px;
                border: 1px solid #e1e1e1;
                font-weight: bold;
            }
            QTableWidget::item {
                color: #111111;
                padding: 4px;
            }
            QTableWidget::item:selected {
                background-color: #111111;
                color: #ffffff;
            }
        """)

        layout.addWidget(self.pigment_table)

        # Add default pigments
        self.add_default_pigments()

        # Buttons
        btn_layout = QHBoxLayout()

        self.add_pigment_btn = QPushButton("+ Add Pigment")
        self.add_pigment_btn.setStyleSheet("""
            QPushButton {
                background-color: #111111;
                color: #ffffff;
                border: 1px solid #111111;
                padding: 6px 12px;
                border-radius: 2px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #000000;
            }
        """)
        self.add_pigment_btn.clicked.connect(self.add_pigment_row)
        btn_layout.addWidget(self.add_pigment_btn)

        self.edit_pigment_btn = QPushButton("Edit Selected")
        self.edit_pigment_btn.setStyleSheet("""
            QPushButton {
                background-color: #ffffff;
                color: #111111;
                border: 1px solid #111111;
                padding: 6px 12px;
                border-radius: 2px;
            }
            QPushButton:hover {
                background-color: #efefef;
            }
        """)
        btn_layout.addWidget(self.edit_pigment_btn)

        self.remove_pigment_btn = QPushButton("Remove")
        self.remove_pigment_btn.setStyleSheet("""
            QPushButton {
                background-color: #ffffff;
                color: #111111;
                border: 1px solid #111111;
                padding: 6px 12px;
                border-radius: 2px;
            }
            QPushButton:hover {
                background-color: #efefef;
            }
        """)
        self.remove_pigment_btn.clicked.connect(self.remove_pigment_row)
        btn_layout.addWidget(self.remove_pigment_btn)

        btn_layout.addStretch()

        layout.addLayout(btn_layout)

        return group

    def add_default_pigments(self):
        """Add default pigment configurations"""
        default_pigments = [
            ("CLA", "Chl-a", "CLA_IPPC", "CLA", "14900"),
            ("CHL", "Chl-b", "CHL_IPPC", "CHL", "15674"),
        ]

        for residue, pigment, tresp, cdc, evac in default_pigments:
            row = self.pigment_table.rowCount()
            self.pigment_table.insertRow(row)
            self.pigment_table.setItem(row, 0, QTableWidgetItem(residue))
            self.pigment_table.setItem(row, 1, QTableWidgetItem(pigment))
            self.pigment_table.setItem(row, 2, QTableWidgetItem(tresp))
            self.pigment_table.setItem(row, 3, QTableWidgetItem(cdc))
            self.pigment_table.setItem(row, 4, QTableWidgetItem(evac))

    def add_pigment_row(self):
        """Add a new pigment row"""
        row = self.pigment_table.rowCount()
        self.pigment_table.insertRow(row)
        self.pigment_table.setItem(row, 0, QTableWidgetItem(""))
        self.pigment_table.setItem(row, 1, QTableWidgetItem(""))
        self.pigment_table.setItem(row, 2, QTableWidgetItem(""))
        self.pigment_table.setItem(row, 3, QTableWidgetItem(""))
        self.pigment_table.setItem(row, 4, QTableWidgetItem(""))

    def remove_pigment_row(self):
        """Remove selected pigment row"""
        current_row = self.pigment_table.currentRow()
        if current_row >= 0:
            self.pigment_table.removeRow(current_row)
            self.pigment_removed.emit(current_row)

    def on_parameters_changed(self):
        """Emit signal when parameters change"""
        params = self.get_parameters()
        self.parameters_changed.emit(params)

    def get_parameters(self) -> Dict[str, Any]:
        """Get current parameter values"""
        return {
            'coupling_method': 'tresp' if self.tresp_radio.isChecked() else 'dipole',
            'dielectric_cdc': self.dielectric_cdc_spin.value(),
            'dielectric_tresp': self.dielectric_tresp_spin.value(),
            'oscillator_strength': self.osc_strength_spin.value(),
            'temperature': self.temperature_spin.value(),
            'pigments': self.get_pigment_configs()
        }

    def get_pigment_configs(self):
        """Get pigment configurations from table"""
        configs = []
        for row in range(self.pigment_table.rowCount()):
            residue = self.pigment_table.item(row, 0)
            pigment = self.pigment_table.item(row, 1)
            tresp = self.pigment_table.item(row, 2)
            cdc = self.pigment_table.item(row, 3)
            evac = self.pigment_table.item(row, 4)

            if residue and pigment:
                config = {
                    'residue': residue.text(),
                    'pigment': pigment.text(),
                    'tresp_dict': tresp.text() if tresp else '',
                    'cdc_dict': cdc.text() if cdc else '',
                    'e_vac': int(evac.text()) if evac and evac.text() else 0
                }
                configs.append(config)

        return configs

    def set_parameters(self, params: Dict[str, Any]):
        """Set parameter values"""
        if 'coupling_method' in params:
            if params['coupling_method'] == 'tresp':
                self.tresp_radio.setChecked(True)
            else:
                self.dipole_radio.setChecked(True)

        if 'dielectric_cdc' in params:
            self.dielectric_cdc_spin.setValue(params['dielectric_cdc'])

        if 'dielectric_tresp' in params:
            self.dielectric_tresp_spin.setValue(params['dielectric_tresp'])

        if 'oscillator_strength' in params:
            self.osc_strength_spin.setValue(params['oscillator_strength'])

        if 'temperature' in params:
            self.temperature_spin.setValue(params['temperature'])
