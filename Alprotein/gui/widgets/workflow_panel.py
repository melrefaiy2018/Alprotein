"""
Workflow Panel Widget for Alprotein GUI

This widget provides a step-by-step workflow interface guiding users through:
1. Structure Input (PDB file loading)
2. Pigment Setup (configuration)
3. Site Energies (calculation)
4. Hamiltonian (construction & diagonalization)
5. Spectra (absorption calculation)
"""

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QFileDialog, QGroupBox, QCheckBox, QFrame, QScrollArea,
    QSizePolicy, QSpacerItem
)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont, QIcon
import sys


class WorkflowStepWidget(QWidget):
    """Individual workflow step widget"""

    def __init__(self, step_number, title, parent=None):
        super().__init__(parent)
        self.step_number = step_number
        self.title = title
        self.is_complete = False
        self.is_active = False
        self.setup_ui()

    def setup_ui(self):
        """Setup the UI for this workflow step"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(10)

        # Header with step number and title
        header_layout = QHBoxLayout()

        # Step indicator (circle with number)
        self.step_indicator = QLabel(f"{self.step_number}")
        self.step_indicator.setFixedSize(28, 28)
        self.step_indicator.setAlignment(Qt.AlignCenter)

        header_layout.addWidget(self.step_indicator)

        # Title
        self.title_label = QLabel(self.title)
        title_font = QFont()
        title_font.setBold(True)
        title_font.setPointSize(11)
        self.title_label.setFont(title_font)
        header_layout.addWidget(self.title_label)

        # Status indicator
        self.status_label = QLabel("")
        self.status_label.setStyleSheet("color: #1f7a44; font-size: 16px;")
        header_layout.addWidget(self.status_label)

        header_layout.addStretch()

        layout.addLayout(header_layout)

        # Content area (will be populated by subclasses)
        self.content_widget = QWidget()
        self.content_layout = QVBoxLayout(self.content_widget)
        self.content_layout.setContentsMargins(32, 0, 0, 0)  # Indent content
        layout.addWidget(self.content_widget)

        # Separator line
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setStyleSheet("background-color: #e1e1e1;")
        separator.setMaximumHeight(1)
        layout.addWidget(separator)

        # Update style now that all widgets are created
        self.update_style()

    def update_style(self):
        """Update the visual style based on state"""
        if self.is_complete:
            bg_color = "#111111"
            text_color = "#ffffff"
            self.status_label.setText("✓")
        elif self.is_active:
            bg_color = "#3a3a3a"
            text_color = "#ffffff"
            self.status_label.setText("")
        else:
            bg_color = "#e1e1e1"
            text_color = "#6b6b6b"
            self.status_label.setText("")

        self.step_indicator.setStyleSheet(f"""
            QLabel {{
                background-color: {bg_color};
                color: {text_color};
                border-radius: 2px;
                font-weight: bold;
                font-size: 12px;
            }}
        """)

    def set_complete(self, complete=True):
        """Mark this step as complete"""
        self.is_complete = complete
        self.update_style()

    def set_active(self, active=True):
        """Mark this step as active"""
        self.is_active = active
        self.update_style()


class StructureInputStep(WorkflowStepWidget):
    """Step 1: Structure Input"""

    file_selected = pyqtSignal(str)  # file_path
    file_cleared = pyqtSignal()
    charges_selected = pyqtSignal(str)  # charges_path

    def __init__(self, parent=None):
        super().__init__(1, "Structure Input", parent)
        self.setup_content()

    def setup_content(self):
        """Setup step-specific content"""
        # PDB file section
        pdb_label = QLabel("Load PDB File")
        pdb_label.setStyleSheet("color: #111111; background-color: transparent; font-weight: bold;")
        self.content_layout.addWidget(pdb_label)

        self.file_path_label = QLabel("No file loaded")
        self.file_path_label.setStyleSheet("color: #6b6b6b; background-color: transparent; font-style: italic;")
        self.content_layout.addWidget(self.file_path_label)

        # Browse button
        btn_layout = QHBoxLayout()
        self.browse_btn = QPushButton("📁 Browse...")
        self.browse_btn.clicked.connect(self.browse_file)
        self.browse_btn.setStyleSheet("""
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
        btn_layout.addWidget(self.browse_btn)

        self.clear_btn = QPushButton("Clear")
        self.clear_btn.clicked.connect(self.clear_file)
        self.clear_btn.setEnabled(False)
        self.clear_btn.setStyleSheet("""
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
            QPushButton:disabled {
                background-color: #f3f3f3;
                color: #b0b0b0;
                border: 1px solid #e1e1e1;
            }
        """)
        btn_layout.addWidget(self.clear_btn)
        btn_layout.addStretch()

        self.content_layout.addLayout(btn_layout)

        # Optional charges section
        self.content_layout.addSpacing(10)
        charges_label = QLabel("Load Charges (Optional PQR)")
        charges_label.setStyleSheet("color: #111111; background-color: transparent;")
        self.content_layout.addWidget(charges_label)

        self.charges_path_label = QLabel("No charges loaded")
        self.charges_path_label.setStyleSheet("color: #6b6b6b; font-style: italic; font-size: 10px;")
        self.content_layout.addWidget(self.charges_path_label)

        self.browse_charges_btn = QPushButton("Browse...")
        self.browse_charges_btn.clicked.connect(self.browse_charges)
        self.browse_charges_btn.setEnabled(False)
        self.browse_charges_btn.setStyleSheet("""
            QPushButton {
                background-color: #ffffff;
                color: #111111;
                border: 1px solid #111111;
                padding: 4px 10px;
                border-radius: 2px;
                font-size: 11px;
            }
            QPushButton:hover {
                background-color: #efefef;
            }
            QPushButton:disabled {
                background-color: #f3f3f3;
                color: #b0b0b0;
                border: 1px solid #e1e1e1;
            }
        """)
        self.content_layout.addWidget(self.browse_charges_btn)

    def browse_file(self):
        """Open file browser for PDB file"""
        options = QFileDialog.DontUseNativeDialog if sys.platform == 'darwin' else QFileDialog.Options()

        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Load PDB File",
            "",
            "PDB Files (*.pdb);;All Files (*)",
            options=options
        )

        if file_path:
            self.set_file(file_path)
            self.file_selected.emit(file_path)

    def set_file(self, file_path):
        """Set the loaded file path"""
        import os
        self.file_path_label.setText(os.path.basename(file_path))
        self.file_path_label.setStyleSheet("color: #111111; background-color: transparent;")
        self.clear_btn.setEnabled(True)
        self.browse_charges_btn.setEnabled(True)
        self.set_complete(True)
        self.set_active(False)

    def clear_file(self):
        """Clear the loaded file"""
        self.file_path_label.setText("No file loaded")
        self.file_path_label.setStyleSheet("color: #6b6b6b; font-style: italic;")
        self.clear_btn.setEnabled(False)
        self.browse_charges_btn.setEnabled(False)
        self.set_complete(False)
        self.file_cleared.emit()

    def browse_charges(self):
        """Browse for PQR charges file"""
        options = QFileDialog.DontUseNativeDialog if sys.platform == 'darwin' else QFileDialog.Options()

        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Load PQR Charges File",
            "",
            "PQR Files (*.pqr);;All Files (*)",
            options=options
        )

        if file_path:
            import os
            self.charges_path_label.setText(os.path.basename(file_path))
            self.charges_path_label.setStyleSheet(
                "color: #111111; background-color: transparent; font-size: 10px;"
            )
            self.charges_selected.emit(file_path)


class PigmentSetupStep(WorkflowStepWidget):
    """Step 2: Pigment Setup"""

    configure_requested = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(2, "Pigment Setup", parent)
        self.setup_content()

    def setup_content(self):
        """Setup step-specific content"""
        self.info_label = QLabel("Configure pigment properties")
        self.info_label.setStyleSheet("color: #6b6b6b;")
        self.content_layout.addWidget(self.info_label)

        self.pigment_count_label = QLabel("Pigments: Not configured")
        self.pigment_count_label.setStyleSheet("color: #6b6b6b; font-size: 10px;")
        self.content_layout.addWidget(self.pigment_count_label)

        self.configure_btn = QPushButton("⚙ Configure...")
        self.configure_btn.clicked.connect(self.configure_requested.emit)
        self.configure_btn.setEnabled(False)
        self.configure_btn.setStyleSheet("""
            QPushButton {
                background-color: #111111;
                color: #ffffff;
                border: 1px solid #111111;
                padding: 6px 12px;
                border-radius: 2px;
            }
            QPushButton:hover {
                background-color: #000000;
            }
            QPushButton:disabled {
                background-color: #f3f3f3;
                color: #b0b0b0;
                border: 1px solid #e1e1e1;
            }
        """)
        self.content_layout.addWidget(self.configure_btn)

    def enable(self, enabled=True):
        """Enable or disable this step"""
        self.configure_btn.setEnabled(enabled)
        if enabled:
            self.set_active(True)

    def set_pigment_info(self, count, pigment_types):
        """Update pigment information"""
        if count > 0:
            types_str = ", ".join(set(pigment_types.values()))
            self.pigment_count_label.setText(f"Pigments: {count} ({types_str})")
            self.pigment_count_label.setStyleSheet(
                "color: #111111; background-color: transparent; font-size: 10px;"
            )
            self.set_complete(True)
            self.set_active(False)


class CalculationStep(WorkflowStepWidget):
    """Generic calculation step"""

    calculate_requested = pyqtSignal()

    def __init__(self, step_number, title, description, parent=None):
        self.description = description
        super().__init__(step_number, title, parent)
        self.setup_content()

    def setup_content(self):
        """Setup step-specific content"""
        desc_label = QLabel(self.description)
        desc_label.setStyleSheet("color: #6b6b6b;")
        desc_label.setWordWrap(True)
        self.content_layout.addWidget(desc_label)

        self.calculate_btn = QPushButton(f"▶ Calculate")
        self.calculate_btn.clicked.connect(self.calculate_requested.emit)
        self.calculate_btn.setEnabled(False)
        self.calculate_btn.setStyleSheet("""
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
            QPushButton:disabled {
                background-color: #f3f3f3;
                color: #b0b0b0;
                border: 1px solid #e1e1e1;
            }
        """)
        self.content_layout.addWidget(self.calculate_btn)

    def enable(self, enabled=True):
        """Enable or disable this step"""
        self.calculate_btn.setEnabled(enabled)
        if enabled:
            self.set_active(True)


class WorkflowPanel(QWidget):
    """
    Left-side workflow panel that guides users through the calculation process
    """

    # Signals
    file_selected = pyqtSignal(str)
    charges_selected = pyqtSignal(str)
    configure_pigments = pyqtSignal()
    calculate_site_energies = pyqtSignal()
    construct_hamiltonian = pyqtSignal()
    diagonalize_hamiltonian = pyqtSignal()
    calculate_spectrum = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setup_ui()

    def setup_ui(self):
        """Setup the UI"""
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)

        # Header
        header = QLabel("WORKFLOW PANEL")
        header.setAlignment(Qt.AlignCenter)
        header.setStyleSheet("""
            QLabel {
                background-color: #f3f3f3;
                color: #111111;
                padding: 14px;
                font-size: 12px;
                font-weight: bold;
                letter-spacing: 1px;
                border-bottom: 1px solid #e1e1e1;
            }
        """)
        main_layout.addWidget(header)

        # Scroll area for steps
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setStyleSheet("""
            QScrollArea {
                border: none;
                background-color: #ffffff;
            }
        """)

        # Container widget for steps
        container = QWidget()
        self.steps_layout = QVBoxLayout(container)
        self.steps_layout.setContentsMargins(16, 16, 16, 16)
        self.steps_layout.setSpacing(0)

        # Create workflow steps
        self.step1 = StructureInputStep()
        self.step1.file_selected.connect(self.file_selected.emit)
        self.step1.file_selected.connect(self.on_file_loaded)
        self.step1.file_cleared.connect(self.on_file_cleared)
        self.step1.charges_selected.connect(self.charges_selected.emit)
        self.steps_layout.addWidget(self.step1)

        self.step2 = PigmentSetupStep()
        self.step2.configure_requested.connect(self.configure_pigments.emit)
        self.steps_layout.addWidget(self.step2)

        self.step3 = CalculationStep(3, "Site Energies", "Calculate CDC site energies")
        self.step3.calculate_requested.connect(self.calculate_site_energies.emit)
        self.steps_layout.addWidget(self.step3)

        self.step4_construct = CalculationStep(4, "Hamiltonian", "Construct Hamiltonian matrix")
        self.step4_construct.calculate_requested.connect(self.construct_hamiltonian.emit)
        self.steps_layout.addWidget(self.step4_construct)

        # Add diagonalize button to Hamiltonian step
        self.diag_btn = QPushButton("🔢 Diagonalize")
        self.diag_btn.clicked.connect(self.diagonalize_hamiltonian.emit)
        self.diag_btn.setEnabled(False)
        self.diag_btn.setStyleSheet("""
            QPushButton {
                background-color: #ffffff;
                color: #111111;
                border: 1px solid #111111;
                padding: 6px 12px;
                border-radius: 2px;
                margin-left: 32px;
            }
            QPushButton:hover {
                background-color: #efefef;
            }
            QPushButton:disabled {
                background-color: #f3f3f3;
                color: #b0b0b0;
                border: 1px solid #e1e1e1;
            }
        """)
        self.step4_construct.content_layout.addWidget(self.diag_btn)

        self.step5 = CalculationStep(5, "Spectra", "Calculate absorption spectrum")
        self.step5.calculate_requested.connect(self.calculate_spectrum.emit)
        self.steps_layout.addWidget(self.step5)

        # Add stretch to push steps to top
        self.steps_layout.addStretch()

        scroll.setWidget(container)
        main_layout.addWidget(scroll)

        # Set initial state
        self.step1.set_active(True)

    def on_file_loaded(self, file_path):
        """Handle file loaded"""
        self.step2.enable(True)

    def on_file_cleared(self):
        """Handle file cleared"""
        self.reset_workflow()

    def reset_workflow(self):
        """Reset workflow to initial state"""
        self.step1.set_active(True)
        self.step2.enable(False)
        self.step2.set_complete(False)
        self.step2.set_active(False)
        self.step3.enable(False)
        self.step3.set_complete(False)
        self.step3.set_active(False)
        self.step4_construct.enable(False)
        self.step4_construct.set_complete(False)
        self.step4_construct.set_active(False)
        self.diag_btn.setEnabled(False)
        self.step5.enable(False)
        self.step5.set_complete(False)
        self.step5.set_active(False)

    def set_pigment_info(self, count, pigment_types):
        """Update pigment information"""
        self.step2.set_pigment_info(count, pigment_types)
        self.step3.enable(True)

    def set_site_energies_complete(self):
        """Mark site energies as complete"""
        self.step3.set_complete(True)
        self.step3.set_active(False)
        self.step4_construct.enable(True)

    def set_hamiltonian_constructed(self):
        """Mark Hamiltonian construction as complete"""
        self.step4_construct.set_complete(True)
        self.step4_construct.set_active(False)
        self.diag_btn.setEnabled(True)

    def set_hamiltonian_diagonalized(self):
        """Mark Hamiltonian diagonalization as complete"""
        self.step5.enable(True)

    def set_spectrum_complete(self):
        """Mark spectrum calculation as complete"""
        self.step5.set_complete(True)
        self.step5.set_active(False)
