"""
Calculation Control Panel Widget
"""

from typing import Dict, Any
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton, 
                            QLabel, QGroupBox, QSpinBox, QDoubleSpinBox, 
                            QCheckBox, QComboBox, QTextEdit, QFrame, QSlider,
                            QFormLayout, QSizePolicy)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont


class CalculationPanel(QWidget):
    """
    Panel for controlling site energy calculations and CDC analysis
    """
    settings_changed = pyqtSignal(dict)  # settings dict
    
    def __init__(self):
        super().__init__()
        self.setup_ui()
        self.connect_signals()
    
    def setup_ui(self):
        """Setup the UI components"""
        layout = QVBoxLayout(self)
        layout.setSpacing(16)
        
        # System Status Group
        # self.create_system_status_group(layout)
        
        # Calculation Settings Group
        self.create_calculation_settings_group(layout)
        
        # Actions Group
        self.create_actions_group(layout)
        
        layout.addStretch()
    
    def create_system_status_group(self, parent_layout):
        """Create system status display group"""
        group_box = QGroupBox("🔍 System Status")
        group_box.setProperty("class", "card")
        layout = QVBoxLayout(group_box)
        layout.setContentsMargins(16, 24, 16, 16)
        
        # Status text area
        self.status_text = QTextEdit()
        self.status_text.setMaximumHeight(120)
        self.status_text.setReadOnly(True)
        # Use global style for QTextEdit if available, or keep minimal styling
        self.status_text.setStyleSheet("""
            QTextEdit {
                background-color: #f9fafb;
                border: 1px solid #e5e7eb;
                border-radius: 6px;
                padding: 8px;
                font-family: "Courier New", monospace;
                font-size: 12px;
            }
        """)
        self.status_text.setText("No structure loaded")
        layout.addWidget(self.status_text)
        
        parent_layout.addWidget(group_box)
    
    def create_calculation_settings_group(self, parent_layout):
        """Create calculation settings group"""
        group_box = QGroupBox("⚙️ Calculation Settings")
        group_box.setProperty("class", "card")
        layout = QFormLayout(group_box)
        layout.setContentsMargins(16, 24, 16, 16)
        layout.setSpacing(12)
        
        # Dielectric constant
        self.dielectric_spinbox = QDoubleSpinBox()
        self.dielectric_spinbox.setRange(1.0, 50000.0)  # Expanded range for high values
        self.dielectric_spinbox.setSingleStep(0.1)
        self.dielectric_spinbox.setValue(2.5)
        self.dielectric_spinbox.setDecimals(2)  # Allow up to 2 decimal places
        self.dielectric_spinbox.setSuffix("")
        self.dielectric_spinbox.setMinimumWidth(150)  # Make it wider
        layout.addRow("Dielectric Constant:", self.dielectric_spinbox)
        
        # E_0a
        self.e0a_spinbox = QDoubleSpinBox()
        self.e0a_spinbox.setRange(-999999.0, 999999.0)  # Much wider range
        self.e0a_spinbox.setSingleStep(0.1)
        self.e0a_spinbox.setValue(14900.0)  # Default value
        self.e0a_spinbox.setDecimals(2)  # Allow up to 2 decimal places
        self.e0a_spinbox.setSuffix(" cm⁻¹")
        self.e0a_spinbox.setMinimumWidth(150)  # Make it wider
        layout.addRow("E₀a (cm⁻¹):", self.e0a_spinbox)
        
        # E_0b
        self.e0b_spinbox = QDoubleSpinBox()
        self.e0b_spinbox.setRange(-999999.0, 999999.0)  # Much wider range
        self.e0b_spinbox.setSingleStep(0.1)
        self.e0b_spinbox.setValue(15674.0)  # Default value
        self.e0b_spinbox.setDecimals(2)  # Allow up to 2 decimal places
        self.e0b_spinbox.setSuffix(" cm⁻¹")
        self.e0b_spinbox.setMinimumWidth(150)  # Make it wider
        layout.addRow("E₀b (cm⁻¹):", self.e0b_spinbox)
        
        # Method selection
        self.method_combo = QComboBox()
        self.method_combo.addItems(["CDC (Charge Density Coupling)"])
        self.method_combo.setCurrentText("CDC (Charge Density Coupling)")
        self.method_combo.setEnabled(False)  # Only CDC supported for now
        layout.addRow("Calculation Method:", self.method_combo)
        
        # CDC cutoff for analysis
        self.cdc_cutoff_spinbox = QDoubleSpinBox()
        self.cdc_cutoff_spinbox.setRange(0.0, 999999.0)  # Much wider range
        self.cdc_cutoff_spinbox.setSingleStep(1.0)
        self.cdc_cutoff_spinbox.setValue(20.0)
        self.cdc_cutoff_spinbox.setDecimals(2)  # Allow up to 2 decimal places
        self.cdc_cutoff_spinbox.setSuffix(" cm⁻¹")
        self.cdc_cutoff_spinbox.setMinimumWidth(150)  # Make it wider
        layout.addRow("CDC Analysis Cutoff:", self.cdc_cutoff_spinbox)
        
        # Advanced options
        self.create_advanced_options(layout)
        
        parent_layout.addWidget(group_box)
    
    def create_advanced_options(self, parent_layout):
        """Create advanced calculation options"""
        # Separator
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        parent_layout.addRow(separator)
        
        # Advanced label
        advanced_label = QLabel("Advanced Options:")
        advanced_label.setProperty("class", "text-primary")
        font = advanced_label.font()
        font.setBold(True)
        advanced_label.setFont(font)
        parent_layout.addRow(advanced_label)
        
        # Include water molecules
        self.include_water_checkbox = QCheckBox("Include water molecules")
        self.include_water_checkbox.setChecked(True)
        parent_layout.addRow("", self.include_water_checkbox)
        
        # Generate plots
        self.generate_plots_checkbox = QCheckBox("Generate visualization plots")
        self.generate_plots_checkbox.setChecked(True)
        parent_layout.addRow("", self.generate_plots_checkbox)
        
        # Detailed output
        self.detailed_output_checkbox = QCheckBox("Generate detailed reports")
        self.detailed_output_checkbox.setChecked(True)
        parent_layout.addRow("", self.detailed_output_checkbox)
    
    def create_actions_group(self, parent_layout):
        """Create action buttons group"""
        # Actions moved to separate tab in control sidebar
        pass
    
    def connect_signals(self):
        """Connect widget signals"""
        # Settings changes
        self.dielectric_spinbox.valueChanged.connect(self.emit_settings_changed)
        self.e0a_spinbox.valueChanged.connect(self.emit_settings_changed)
        self.e0b_spinbox.valueChanged.connect(self.emit_settings_changed)
        self.cdc_cutoff_spinbox.valueChanged.connect(self.emit_settings_changed)
        self.include_water_checkbox.toggled.connect(self.emit_settings_changed)
        self.generate_plots_checkbox.toggled.connect(self.emit_settings_changed)
        self.detailed_output_checkbox.toggled.connect(self.emit_settings_changed)
    
    def emit_settings_changed(self):
        """Emit settings changed signal with current settings"""
        settings = self.get_current_settings()
        self.settings_changed.emit(settings)
    
    def get_current_settings(self) -> Dict[str, Any]:
        """Get current calculation settings"""
        return {
            'dielectric_constant': self.dielectric_spinbox.value(),
            'e0a': self.e0a_spinbox.value(),
            'e0b': self.e0b_spinbox.value(),
            'method': self.method_combo.currentText(),
            'cdc_cutoff': self.cdc_cutoff_spinbox.value(),
            'include_water': self.include_water_checkbox.isChecked(),
            'generate_plots': self.generate_plots_checkbox.isChecked(),
            'detailed_output': self.detailed_output_checkbox.isChecked()
        }
    
    def set_system_info(self, validation: Dict[str, Any]):
        """Update system status display"""
        if not validation:
            self.status_text.setText("No structure loaded")
            self.system_loaded = False
            self.update_button_states()
            return
        
        # Format validation info
        status_lines = [
            f"✅ Structure loaded successfully",
            f"📊 Pigments: {validation['total_pigments']}",
            f"🧬 Background atoms: {validation['background_atoms']}",
            "",
            "Pigment Types:"
        ]
        
        for ptype, count in validation['pigment_types'].items():
            status_lines.append(f"  • {ptype}: {count}")
        
        status_lines.append("")
        status_lines.append("Calculation Readiness:")
        
        readiness = validation['calculation_readiness']
        total = validation['total_pigments']
        
        status_lines.extend([
            f"  • Site Energy: {readiness['site_energy_ready']}/{total}",
            f"  • Coupling: {readiness['coupling_ready']}/{total}",
            f"  • Dipoles: {readiness['dipole_ready']}/{total}",
            f"  • Fully Ready: {readiness['fully_ready']}/{total}"
        ])
        
        if validation.get('warnings'):
            status_lines.append("")
            status_lines.append("⚠️ Warnings:")
            for warning in validation['warnings'][:5]:  # Show first 5 warnings
                status_lines.append(f"  • {warning}")
            if len(validation['warnings']) > 5:
                status_lines.append(f"  • ... and {len(validation['warnings']) - 5} more")
        
        self.status_text.setText("\n".join(status_lines))
        # Keep track of system loaded state for external components
        self.system_loaded = True
    
    
