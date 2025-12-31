"""
Tools Panel Widget for Scientific Workbench Layout

This panel provides project information, settings, and results tracking
in a professional sidebar format.
"""

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QScrollArea, QFrame, QLineEdit, QComboBox,
    QDoubleSpinBox, QSpinBox, QRadioButton, QCheckBox,
    QSizePolicy, QButtonGroup, QProgressBar
)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont, QCursor
from typing import Dict, Any


class DragDropArea(QLabel):
    """Drag-and-drop area for PDB files"""
    file_dropped = pyqtSignal(str)

    def __init__(self):
        super().__init__("🧬 Drag PDB file here\nor click to browse")
        self.setAcceptDrops(True)
        self.setAlignment(Qt.AlignCenter)
        self.setMinimumHeight(80)
        self.setStyleSheet("""
            QLabel {
                border: 1px dashed #e1e1e1;
                border-radius: 2px;
                background-color: #f7f7f7;
                color: #6b6b6b;
                font-size: 13px;
                padding: 10px;
            }
            QLabel:hover {
                border-color: #111111;
                background-color: #efefef;
                color: #111111;
            }
        """)
        self.setCursor(QCursor(Qt.PointingHandCursor))

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()
            self.setStyleSheet("""
                QLabel {
                    border: 2px solid #1f7a44;
                    border-radius: 2px;
                    background-color: #eaf3ed;
                    color: #1f7a44;
                }
            """)

    def dragLeaveEvent(self, event):
        self.setStyleSheet("""
            QLabel {
                border: 1px dashed #e1e1e1;
                border-radius: 2px;
                background-color: #f7f7f7;
                color: #6b6b6b;
            }
        """)

    def dropEvent(self, event):
        urls = event.mimeData().urls()
        if urls:
            file_path = urls[0].toLocalFile()
            if file_path.endswith('.pdb'):
                self.file_dropped.emit(file_path)
            self.dragLeaveEvent(event)  # Reset style

    def mousePressEvent(self, event):
        self.file_dropped.emit("")  # Signal to open file dialog


class ToolsPanel(QWidget):
    """
    Right sidebar tools panel for the scientific workbench
    """

    # Signals
    open_project = pyqtSignal()
    show_properties = pyqtSignal()
    parameters_changed = pyqtSignal(dict)
    run_calculation = pyqtSignal(str)  # calculation type
    view_results = pyqtSignal(str)  # result type
    file_dropped = pyqtSignal(str)  # PDB file path from drag-drop

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumWidth(280)
        self.setMaximumWidth(350)
        self.results_items = {}
        self.setup_ui()

    def setup_ui(self):
        """Setup the UI"""
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(16, 16, 16, 16)
        main_layout.setSpacing(16)

        # Create scroll area for all tools
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setFrameShape(QFrame.NoFrame)

        container = QWidget()
        container_layout = QVBoxLayout(container)
        container_layout.setSpacing(18)

        # Project section
        self.project_group = self.create_project_section()
        container_layout.addWidget(self.project_group)

        # Settings section
        self.settings_group = self.create_settings_section()
        container_layout.addWidget(self.settings_group)

        container_layout.addStretch()

        scroll.setWidget(container)
        main_layout.addWidget(scroll)

        # Apply styling with improved QGroupBox title visibility
        self.setStyleSheet("""
            ToolsPanel {
                background-color: #f5f5f5;
            }
            QScrollArea {
                background-color: #f5f5f5;
                border: none;
            }
            QScrollArea > QWidget > QWidget {
                background-color: #f5f5f5;
            }
            QGroupBox {
                font-weight: bold;
                font-size: 12px;
                border: 1px solid #e1e1e1;
                border-radius: 4px;
                margin-top: 16px;
                padding: 12px;
                padding-top: 24px;
                background-color: #ffffff;
                color: #111111;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                subcontrol-position: top left;
                left: 12px;
                top: 4px;
                padding: 2px 8px;
                background-color: #ffffff;
                color: #111111;
                font-size: 12px;
                font-weight: bold;
            }
            QLabel {
                color: #111111;
                background-color: transparent;
            }
            QPushButton {
                background-color: #111111;
                color: #ffffff;
                border: 1px solid #111111;
                padding: 8px 12px;
                border-radius: 2px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #000000;
            }
            QPushButton:pressed {
                background-color: #000000;
            }
            QLineEdit, QComboBox, QDoubleSpinBox, QSpinBox {
                padding: 6px;
                border: 1px solid #e1e1e1;
                border-radius: 2px;
                background-color: #ffffff;
                color: #111111;
            }
            QLineEdit:focus, QComboBox:focus, QDoubleSpinBox:focus, QSpinBox:focus {
                border: 2px solid #111111;
            }
            QCheckBox {
                color: #111111;
                background-color: transparent;
            }
            QRadioButton {
                color: #111111;
                background-color: transparent;
            }
            QProgressBar {
                border: 1px solid #e1e1e1;
                border-radius: 2px;
                text-align: center;
                background-color: #f3f3f3;
                color: #111111;
            }
            QProgressBar::chunk {
                background-color: #111111;
                border-radius: 2px;
            }
        """)

    def create_project_section(self):
        """Create project information section"""
        group = QGroupBox("📁 PROJECT")
        group.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                font-size: 12px;
                border: 1px solid #d0d0d0;
                border-radius: 4px;
                margin-top: 16px;
                padding: 12px;
                padding-top: 24px;
                background-color: #ffffff;
                color: #111111;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                subcontrol-position: top left;
                left: 12px;
                top: 4px;
                padding: 2px 8px;
                background-color: #ffffff;
                color: #111111;
                font-size: 12px;
                font-weight: bold;
            }
        """)
        layout = QVBoxLayout(group)
        layout.setSpacing(10)

        # Drag-drop area
        self.drag_drop_area = DragDropArea()
        self.drag_drop_area.file_dropped.connect(self.file_dropped.emit)
        layout.addWidget(self.drag_drop_area)

        # Separator
        sep = QFrame()
        sep.setFrameShape(QFrame.HLine)
        sep.setStyleSheet("background-color: #e1e1e1;")
        layout.addWidget(sep)

        # Project name
        self.project_name_label = QLabel("No project loaded")
        self.project_name_label.setWordWrap(True)
        self.project_name_label.setStyleSheet("font-weight: bold; font-size: 14px; color: #111111;")
        layout.addWidget(self.project_name_label)

        # Separator
        sep = QFrame()
        sep.setFrameShape(QFrame.HLine)
        sep.setStyleSheet("background-color: #e1e1e1;")
        layout.addWidget(sep)

        # Statistics
        stats_layout = QVBoxLayout()
        stats_layout.setSpacing(5)

        self.pigments_label = QLabel("Pigments: -")
        self.pigments_label.setStyleSheet("color: #6b6b6b;")
        stats_layout.addWidget(self.pigments_label)

        self.atoms_label = QLabel("Atoms: -")
        self.atoms_label.setStyleSheet("color: #6b6b6b;")
        stats_layout.addWidget(self.atoms_label)

        self.validation_label = QLabel("Status: Not validated")
        self.validation_label.setStyleSheet("color: #6b6b6b;")
        stats_layout.addWidget(self.validation_label)

        layout.addLayout(stats_layout)

        # Buttons
        btn_layout = QHBoxLayout()

        self.open_btn = QPushButton("Open")
        self.open_btn.clicked.connect(self.open_project.emit)
        btn_layout.addWidget(self.open_btn)

        self.props_btn = QPushButton("Properties")
        self.props_btn.clicked.connect(self.show_properties.emit)
        self.props_btn.setEnabled(False)
        btn_layout.addWidget(self.props_btn)

        layout.addLayout(btn_layout)

        return group

    def create_settings_section(self):
        """Create calculation settings section"""
        group = QGroupBox("⚙ SETTINGS")
        group.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                font-size: 12px;
                border: 1px solid #d0d0d0;
                border-radius: 4px;
                margin-top: 16px;
                padding: 12px;
                padding-top: 24px;
                background-color: #ffffff;
                color: #111111;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                subcontrol-position: top left;
                left: 12px;
                top: 4px;
                padding: 2px 8px;
                background-color: #ffffff;
                color: #111111;
                font-size: 12px;
                font-weight: bold;
            }
        """)
        layout = QVBoxLayout(group)
        layout.setSpacing(12)

        # Calculation header
        calc_header = QLabel("Calculation")
        calc_header.setStyleSheet("font-weight: bold; color: #111111; font-size: 12px;")
        layout.addWidget(calc_header)

        # Dielectric CDC
        cdc_layout = QHBoxLayout()
        cdc_layout.addWidget(QLabel("ε(CDC)"))
        self.dielectric_cdc = QDoubleSpinBox()
        self.dielectric_cdc.setRange(1.0, 80.0)
        self.dielectric_cdc.setValue(2.0)
        self.dielectric_cdc.setDecimals(2)
        self.dielectric_cdc.setSingleStep(0.1)
        self.dielectric_cdc.valueChanged.connect(self.on_settings_changed)
        cdc_layout.addWidget(self.dielectric_cdc)
        layout.addLayout(cdc_layout)

        # Dielectric TrEsp
        tresp_layout = QHBoxLayout()
        tresp_layout.addWidget(QLabel("ε(TrEsp)"))
        self.dielectric_tresp = QDoubleSpinBox()
        self.dielectric_tresp.setRange(1.0, 80.0)
        self.dielectric_tresp.setValue(1.0)
        self.dielectric_tresp.setDecimals(2)
        self.dielectric_tresp.setSingleStep(0.1)
        self.dielectric_tresp.valueChanged.connect(self.on_settings_changed)
        tresp_layout.addWidget(self.dielectric_tresp)
        layout.addLayout(tresp_layout)

        # Oscillator strength
        fosc_layout = QHBoxLayout()
        fosc_layout.addWidget(QLabel("f_osc"))
        self.f_osc = QDoubleSpinBox()
        self.f_osc.setRange(0.1, 2.0)
        self.f_osc.setValue(0.72)
        self.f_osc.setDecimals(3)
        self.f_osc.setSingleStep(0.01)
        self.f_osc.valueChanged.connect(self.on_settings_changed)
        fosc_layout.addWidget(self.f_osc)
        layout.addLayout(fosc_layout)

        # Separator
        sep = QFrame()
        sep.setFrameShape(QFrame.HLine)
        sep.setStyleSheet("background-color: #e1e1e1; margin: 5px 0;")
        layout.addWidget(sep)

        # Reference Energies header
        ref_header = QLabel("Reference Energies (cm⁻¹)")
        ref_header.setStyleSheet("font-weight: bold; color: #111111; font-size: 12px;")
        layout.addWidget(ref_header)

        # Eg,a (CLA)
        ega_layout = QHBoxLayout()
        ega_layout.addWidget(QLabel("Eg,a (CLA)"))
        self.e0a = QDoubleSpinBox()
        self.e0a.setRange(10000, 20000)
        self.e0a.setValue(14900.0)
        self.e0a.setDecimals(1)
        self.e0a.setSuffix(" cm⁻¹")
        self.e0a.valueChanged.connect(self.on_settings_changed)
        ega_layout.addWidget(self.e0a)
        layout.addLayout(ega_layout)

        # Eg,b (CHL)
        egb_layout = QHBoxLayout()
        egb_layout.addWidget(QLabel("Eg,b (CHL)"))
        self.e0b = QDoubleSpinBox()
        self.e0b.setRange(10000, 20000)
        self.e0b.setValue(15674.0)
        self.e0b.setDecimals(1)
        self.e0b.setSuffix(" cm⁻¹")
        self.e0b.valueChanged.connect(self.on_settings_changed)
        egb_layout.addWidget(self.e0b)
        layout.addLayout(egb_layout)

        # Action buttons
        layout.addSpacing(10)

        return group

    def create_progress_section(self):
        """Create progress tracking section"""
        group = QGroupBox("⏱ PROGRESS")
        group.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                font-size: 12px;
                border: 1px solid #d0d0d0;
                border-radius: 4px;
                margin-top: 16px;
                padding: 12px;
                padding-top: 24px;
                background-color: #ffffff;
                color: #111111;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                subcontrol-position: top left;
                left: 12px;
                top: 4px;
                padding: 2px 8px;
                background-color: #ffffff;
                color: #111111;
                font-size: 12px;
                font-weight: bold;
            }
        """)
        layout = QVBoxLayout(group)
        layout.setSpacing(8)

        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(True)
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
                border-radius: 2px;
            }
        """)
        layout.addWidget(self.progress_bar)

        # Status label
        self.progress_status = QLabel("Ready")
        self.progress_status.setStyleSheet("color: #6b6b6b; font-size: 11px;")
        layout.addWidget(self.progress_status)

        # Timing label
        self.progress_timing = QLabel("Elapsed: 0.0s")
        self.progress_timing.setStyleSheet("color: #6b6b6b; font-size: 10px;")
        layout.addWidget(self.progress_timing)

        return group

    def create_results_section(self):
        """Create results tracking section"""
        group = QGroupBox("📊 RESULTS")
        group.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                font-size: 12px;
                border: 1px solid #d0d0d0;
                border-radius: 4px;
                margin-top: 16px;
                padding: 12px;
                padding-top: 24px;
                background-color: #ffffff;
                color: #111111;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                subcontrol-position: top left;
                left: 12px;
                top: 4px;
                padding: 2px 8px;
                background-color: #ffffff;
                color: #111111;
                font-size: 12px;
                font-weight: bold;
            }
        """)
        layout = QVBoxLayout(group)
        layout.setSpacing(10)

        # Results list
        self.results_items = {}

        # Site Energy
        self.site_energy_item = self.create_result_item(
            "Site Energy", "○ Pending", "site_energies"
        )
        layout.addWidget(self.site_energy_item)

        # Hamiltonian
        self.hamiltonian_item = self.create_result_item(
            "Hamiltonian", "○ Pending", "hamiltonian"
        )
        layout.addWidget(self.hamiltonian_item)

        # Spectrum
        self.spectrum_item = self.create_result_item(
            "Spectrum", "○ Pending", "spectrum"
        )
        layout.addWidget(self.spectrum_item)

        # View All button
        layout.addSpacing(10)
        self.view_all_btn = QPushButton("View All Results")
        self.view_all_btn.setStyleSheet("""
            QPushButton {
                background-color: #ffffff;
                color: #111111;
                border: 1px solid #111111;
                padding: 8px;
            }
            QPushButton:hover {
                background-color: #efefef;
            }
        """)
        self.view_all_btn.clicked.connect(lambda: self.view_results.emit("all"))
        layout.addWidget(self.view_all_btn)

        return group

    def create_result_item(self, name, status, result_type):
        """Create a result item widget"""
        item = QFrame()
        item.setFrameShape(QFrame.StyledPanel)
        item.setStyleSheet("""
            QFrame {
                background-color: #f3f3f3;
                border-radius: 2px;
                padding: 8px;
            }
        """)

        item_layout = QVBoxLayout(item)
        item_layout.setContentsMargins(8, 8, 8, 8)
        item_layout.setSpacing(4)

        name_label = QLabel(name)
        name_label.setStyleSheet("font-weight: bold; color: #111111;")
        item_layout.addWidget(name_label)

        status_label = QLabel(status)
        status_label.setStyleSheet("color: #6b6b6b; font-size: 11px;")
        item_layout.addWidget(status_label)

        # Store reference for later updates
        self.results_items[result_type] = {
            'widget': item,
            'status_label': status_label
        }

        return item

    def on_settings_changed(self):
        """Emit settings change signal"""
        params = {
            'dielectric_cdc': self.dielectric_cdc.value(),
            'dielectric_tresp': self.dielectric_tresp.value(),
            'oscillator_strength': self.f_osc.value(),
            'e0a': self.e0a.value(),
            'e0b': self.e0b.value()
        }
        self.parameters_changed.emit(params)

    def update_project_info(self, name, pigments, atoms, validated):
        """Update project information display"""
        self.project_name_label.setText(name)
        self.pigments_label.setText(f"Pigments: {pigments}")
        self.atoms_label.setText(f"Atoms: {atoms:,}")

        if validated:
            self.validation_label.setText("✓ Validated")
            self.validation_label.setStyleSheet("color: #1f7a44; font-weight: bold;")
        else:
            self.validation_label.setText("✗ Not validated")
            self.validation_label.setStyleSheet("color: #8b2b2b; font-weight: bold;")

        self.props_btn.setEnabled(True)
    def update_result_status(self, result_type, status, details=""):
        """Update result status"""
        if result_type in self.results_items:
            item = self.results_items[result_type]

            if status == "complete":
                item['status_label'].setText(f"✓ Complete\n{details}")
                item['status_label'].setStyleSheet("color: #1f7a44; font-weight: bold; font-size: 11px;")
            elif status == "running":
                item['status_label'].setText(f"▶ Running...\n{details}")
                item['status_label'].setStyleSheet("color: #8a6d1f; font-weight: bold; font-size: 11px;")
            elif status == "error":
                item['status_label'].setText(f"✗ Error\n{details}")
                item['status_label'].setStyleSheet("color: #8b2b2b; font-weight: bold; font-size: 11px;")
            else:
                item['status_label'].setText("○ Pending")
                item['status_label'].setStyleSheet("color: #6b6b6b; font-size: 11px;")

    def get_parameters(self):
        """Get current parameter values"""
        return {
            'dielectric_cdc': self.dielectric_cdc.value(),
            'dielectric_tresp': self.dielectric_tresp.value(),
            'oscillator_strength': self.f_osc.value(),
            'e0a': self.e0a.value(),
            'e0b': self.e0b.value(),
        }
