"""
Settings Dialog for calculation parameters
"""

from typing import Dict, Any
from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QGroupBox,
                            QDoubleSpinBox, QSpinBox, QComboBox, QCheckBox,
                            QLabel, QDialogButtonBox, QFormLayout, QTabWidget,
                            QWidget, QTextEdit, QPushButton, QSlider)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont


class SettingsDialog(QDialog):
    """
    Dialog for configuring calculation settings
    """
    
    def __init__(self, current_settings: Dict[str, Any] = None, parent=None):
        super().__init__(parent)
        self.current_settings = current_settings or {}
        
        self.setWindowTitle("Calculation Settings")
        self.setModal(True)
        self.setMinimumSize(600, 500)
        
        self.setup_ui()
        self.load_settings()
    
    def setup_ui(self):
        """Setup the dialog UI"""
        layout = QVBoxLayout(self)
        
        # Title
        title_label = QLabel("⚙️ Calculation Settings")
        title_label.setStyleSheet("font-size: 16px; font-weight: bold; margin-bottom: 10px;")
        layout.addWidget(title_label)
        
        # Create tab widget
        self.tab_widget = QTabWidget()
        
        # Basic settings tab
        self.create_basic_tab()
        self.tab_widget.addTab(self.basic_tab, "Basic")
        
        # Advanced settings tab
        self.create_advanced_tab()
        self.tab_widget.addTab(self.advanced_tab, "Advanced")
        
        # Visualization settings tab
        self.create_visualization_tab()
        self.tab_widget.addTab(self.visualization_tab, "Visualization")
        
        layout.addWidget(self.tab_widget)
        
        # Buttons
        button_layout = QHBoxLayout()
        
        self.reset_button = QPushButton("Reset to Defaults")
        self.reset_button.clicked.connect(self.reset_to_defaults)
        button_layout.addWidget(self.reset_button)
        
        button_layout.addStretch()
        
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        button_layout.addWidget(button_box)
        
        layout.addLayout(button_layout)
    
    def create_basic_tab(self):
        """Create basic settings tab"""
        self.basic_tab = QWidget()
        layout = QVBoxLayout(self.basic_tab)
        
        # Calculation Method
        method_group = QGroupBox("🧮 Calculation Method")
        method_layout = QFormLayout(method_group)
        
        self.method_combo = QComboBox()
        self.method_combo.addItems(["CDC (Charge Density Coupling)"])
        self.method_combo.setEnabled(False)  # Only CDC supported
        method_layout.addRow("Method:", self.method_combo)
        
        layout.addWidget(method_group)
        
        # Environment Parameters
        env_group = QGroupBox("🌍 Environment Parameters")
        env_layout = QFormLayout(env_group)
        
        # Dielectric constant
        self.dielectric_spinbox = QDoubleSpinBox()
        self.dielectric_spinbox.setRange(1.0, 20.0)
        self.dielectric_spinbox.setSingleStep(0.1)
        self.dielectric_spinbox.setValue(2.5)
        self.dielectric_spinbox.setDecimals(1)
        env_layout.addRow("Dielectric Constant:", self.dielectric_spinbox)
        
        # Temperature (future use)
        self.temperature_spinbox = QDoubleSpinBox()
        self.temperature_spinbox.setRange(200.0, 400.0)
        self.temperature_spinbox.setSingleStep(1.0)
        self.temperature_spinbox.setValue(298.15)
        self.temperature_spinbox.setSuffix(" K")
        self.temperature_spinbox.setEnabled(False)  # Not implemented yet
        env_layout.addRow("Temperature:", self.temperature_spinbox)
        
        layout.addWidget(env_group)
        
        # Convergence Criteria
        conv_group = QGroupBox("🎯 Convergence Criteria")
        conv_layout = QFormLayout(conv_group)
        
        # Energy threshold
        self.energy_threshold_spinbox = QDoubleSpinBox()
        self.energy_threshold_spinbox.setRange(0.001, 10.0)
        self.energy_threshold_spinbox.setSingleStep(0.001)
        self.energy_threshold_spinbox.setValue(0.1)
        self.energy_threshold_spinbox.setDecimals(3)
        self.energy_threshold_spinbox.setSuffix(" cm⁻¹")
        conv_layout.addRow("Energy Threshold:", self.energy_threshold_spinbox)
        
        # Maximum iterations
        self.max_iterations_spinbox = QSpinBox()
        self.max_iterations_spinbox.setRange(10, 10000)
        self.max_iterations_spinbox.setSingleStep(10)
        self.max_iterations_spinbox.setValue(1000)
        self.max_iterations_spinbox.setEnabled(False)  # Not used in current implementation
        conv_layout.addRow("Max Iterations:", self.max_iterations_spinbox)
        
        layout.addWidget(conv_group)
        
        layout.addStretch()
    
    def create_advanced_tab(self):
        """Create advanced settings tab"""
        self.advanced_tab = QWidget()
        layout = QVBoxLayout(self.advanced_tab)
        
        # CDC Analysis Parameters
        cdc_group = QGroupBox("📊 CDC Analysis Parameters")
        cdc_layout = QFormLayout(cdc_group)
        
        # CDC cutoff
        self.cdc_cutoff_spinbox = QDoubleSpinBox()
        self.cdc_cutoff_spinbox.setRange(0.1, 1000.0)
        self.cdc_cutoff_spinbox.setSingleStep(1.0)
        self.cdc_cutoff_spinbox.setValue(20.0)
        self.cdc_cutoff_spinbox.setDecimals(1)
        self.cdc_cutoff_spinbox.setSuffix(" cm⁻¹")
        cdc_layout.addRow("Contribution Cutoff:", self.cdc_cutoff_spinbox)
        
        # Distance cutoff
        self.distance_cutoff_spinbox = QDoubleSpinBox()
        self.distance_cutoff_spinbox.setRange(5.0, 100.0)
        self.distance_cutoff_spinbox.setSingleStep(1.0)
        self.distance_cutoff_spinbox.setValue(20.0)
        self.distance_cutoff_spinbox.setDecimals(1)
        self.distance_cutoff_spinbox.setSuffix(" Å")
        cdc_layout.addRow("Distance Cutoff:", self.distance_cutoff_spinbox)
        
        layout.addWidget(cdc_group)
        
        # Calculation Options
        calc_group = QGroupBox("⚡ Calculation Options")
        calc_layout = QVBoxLayout(calc_group)
        
        self.include_water_checkbox = QCheckBox("Include water molecules in background")
        self.include_water_checkbox.setChecked(True)
        calc_layout.addWidget(self.include_water_checkbox)
        
        self.include_ions_checkbox = QCheckBox("Include ions in background")
        self.include_ions_checkbox.setChecked(True)
        calc_layout.addWidget(self.include_ions_checkbox)
        
        self.use_optimized_calc = QCheckBox("Use optimized calculation paths")
        self.use_optimized_calc.setChecked(True)
        calc_layout.addWidget(self.use_optimized_calc)
        
        self.parallel_processing = QCheckBox("Enable parallel processing")
        self.parallel_processing.setChecked(True)
        self.parallel_processing.setEnabled(False)  # Not implemented yet
        calc_layout.addWidget(self.parallel_processing)
        
        layout.addWidget(calc_group)
        
        # Output Options
        output_group = QGroupBox("📄 Output Options")
        output_layout = QVBoxLayout(output_group)
        
        self.detailed_output = QCheckBox("Generate detailed output files")
        self.detailed_output.setChecked(True)
        output_layout.addWidget(self.detailed_output)
        
        self.save_intermediate = QCheckBox("Save intermediate results")
        self.save_intermediate.setChecked(False)
        output_layout.addWidget(self.save_intermediate)
        
        self.verbose_logging = QCheckBox("Enable verbose logging")
        self.verbose_logging.setChecked(False)
        output_layout.addWidget(self.verbose_logging)
        
        layout.addWidget(output_group)
        
        layout.addStretch()
    
    def create_visualization_tab(self):
        """Create visualization settings tab"""
        self.visualization_tab = QWidget()
        layout = QVBoxLayout(self.visualization_tab)
        
        # 3D Visualization
        viz_group = QGroupBox("🎨 3D Visualization")
        viz_layout = QFormLayout(viz_group)
        
        # Default style
        self.default_protein_style = QComboBox()
        self.default_protein_style.addItems(["Cartoon", "Surface", "Line", "Stick"])
        self.default_protein_style.setCurrentText("Cartoon")
        viz_layout.addRow("Default Protein Style:", self.default_protein_style)
        
        self.default_pigment_style = QComboBox()
        self.default_pigment_style.addItems(["Stick", "Ball & Stick", "Sphere", "Line"])
        self.default_pigment_style.setCurrentText("Stick")
        viz_layout.addRow("Default Pigment Style:", self.default_pigment_style)
        
        # Color scheme
        self.default_color_scheme = QComboBox()
        self.default_color_scheme.addItems(["By Element", "By Residue", "By Chain", "By Energy"])
        self.default_color_scheme.setCurrentText("By Element")
        viz_layout.addRow("Default Color Scheme:", self.default_color_scheme)
        
        layout.addWidget(viz_group)
        
        # Plot Settings
        plot_group = QGroupBox("📈 Plot Settings")
        plot_layout = QFormLayout(plot_group)
        
        # Figure size
        self.figure_width_spinbox = QSpinBox()
        self.figure_width_spinbox.setRange(6, 20)
        self.figure_width_spinbox.setValue(12)
        self.figure_width_spinbox.setSuffix(" inches")
        plot_layout.addRow("Figure Width:", self.figure_width_spinbox)
        
        self.figure_height_spinbox = QSpinBox()
        self.figure_height_spinbox.setRange(4, 16)
        self.figure_height_spinbox.setValue(8)
        self.figure_height_spinbox.setSuffix(" inches")
        plot_layout.addRow("Figure Height:", self.figure_height_spinbox)
        
        # DPI
        self.dpi_spinbox = QSpinBox()
        self.dpi_spinbox.setRange(72, 600)
        self.dpi_spinbox.setSingleStep(72)
        self.dpi_spinbox.setValue(300)
        self.dpi_spinbox.setSuffix(" DPI")
        plot_layout.addRow("Resolution:", self.dpi_spinbox)
        
        # Font size
        self.font_size_spinbox = QSpinBox()
        self.font_size_spinbox.setRange(6, 20)
        self.font_size_spinbox.setValue(12)
        self.font_size_spinbox.setSuffix(" pt")
        plot_layout.addRow("Font Size:", self.font_size_spinbox)
        
        layout.addWidget(plot_group)
        
        # Color Settings
        color_group = QGroupBox("🎨 Color Settings")
        color_layout = QVBoxLayout(color_group)
        
        self.use_colorblind_friendly = QCheckBox("Use colorblind-friendly palette")
        self.use_colorblind_friendly.setChecked(False)
        color_layout.addWidget(self.use_colorblind_friendly)
        
        self.high_contrast_mode = QCheckBox("High contrast mode")
        self.high_contrast_mode.setChecked(False)
        color_layout.addWidget(self.high_contrast_mode)
        
        layout.addWidget(color_group)
        
        layout.addStretch()
    
    def load_settings(self):
        """Load current settings into the dialog"""
        if not self.current_settings:
            return
        
        # Basic settings
        if 'dielectric_constant' in self.current_settings:
            self.dielectric_spinbox.setValue(self.current_settings['dielectric_constant'])
        
        if 'energy_threshold' in self.current_settings:
            self.energy_threshold_spinbox.setValue(self.current_settings['energy_threshold'])
        
        # Advanced settings
        if 'cdc_cutoff' in self.current_settings:
            self.cdc_cutoff_spinbox.setValue(self.current_settings['cdc_cutoff'])
        
        if 'distance_cutoff' in self.current_settings:
            self.distance_cutoff_spinbox.setValue(self.current_settings['distance_cutoff'])
        
        if 'include_water' in self.current_settings:
            self.include_water_checkbox.setChecked(self.current_settings['include_water'])
        
        if 'include_ions' in self.current_settings:
            self.include_ions_checkbox.setChecked(self.current_settings['include_ions'])
        
        if 'detailed_output' in self.current_settings:
            self.detailed_output.setChecked(self.current_settings['detailed_output'])
        
        if 'verbose_logging' in self.current_settings:
            self.verbose_logging.setChecked(self.current_settings['verbose_logging'])
        
        # Visualization settings
        if 'default_protein_style' in self.current_settings:
            self.default_protein_style.setCurrentText(self.current_settings['default_protein_style'])
        
        if 'default_pigment_style' in self.current_settings:
            self.default_pigment_style.setCurrentText(self.current_settings['default_pigment_style'])
        
        if 'figure_width' in self.current_settings:
            self.figure_width_spinbox.setValue(self.current_settings['figure_width'])
        
        if 'figure_height' in self.current_settings:
            self.figure_height_spinbox.setValue(self.current_settings['figure_height'])
        
        if 'dpi' in self.current_settings:
            self.dpi_spinbox.setValue(self.current_settings['dpi'])
        
        if 'font_size' in self.current_settings:
            self.font_size_spinbox.setValue(self.current_settings['font_size'])
    
    def reset_to_defaults(self):
        """Reset all settings to defaults"""
        # Basic settings
        self.dielectric_spinbox.setValue(2.5)
        self.temperature_spinbox.setValue(298.15)
        self.energy_threshold_spinbox.setValue(0.1)
        self.max_iterations_spinbox.setValue(1000)
        
        # Advanced settings
        self.cdc_cutoff_spinbox.setValue(20.0)
        self.distance_cutoff_spinbox.setValue(20.0)
        self.include_water_checkbox.setChecked(True)
        self.include_ions_checkbox.setChecked(True)
        self.use_optimized_calc.setChecked(True)
        self.detailed_output.setChecked(True)
        self.save_intermediate.setChecked(False)
        self.verbose_logging.setChecked(False)
        
        # Visualization settings
        self.default_protein_style.setCurrentText("Cartoon")
        self.default_pigment_style.setCurrentText("Stick")
        self.default_color_scheme.setCurrentText("By Element")
        self.figure_width_spinbox.setValue(12)
        self.figure_height_spinbox.setValue(8)
        self.dpi_spinbox.setValue(300)
        self.font_size_spinbox.setValue(12)
        self.use_colorblind_friendly.setChecked(False)
        self.high_contrast_mode.setChecked(False)
    
    def get_settings(self) -> Dict[str, Any]:
        """Get all current settings"""
        return {
            # Basic settings
            'method': self.method_combo.currentText(),
            'dielectric_constant': self.dielectric_spinbox.value(),
            'temperature': self.temperature_spinbox.value(),
            'energy_threshold': self.energy_threshold_spinbox.value(),
            'max_iterations': self.max_iterations_spinbox.value(),
            
            # Advanced settings
            'cdc_cutoff': self.cdc_cutoff_spinbox.value(),
            'distance_cutoff': self.distance_cutoff_spinbox.value(),
            'include_water': self.include_water_checkbox.isChecked(),
            'include_ions': self.include_ions_checkbox.isChecked(),
            'use_optimized_calc': self.use_optimized_calc.isChecked(),
            'parallel_processing': self.parallel_processing.isChecked(),
            'detailed_output': self.detailed_output.isChecked(),
            'save_intermediate': self.save_intermediate.isChecked(),
            'verbose_logging': self.verbose_logging.isChecked(),
            
            # Visualization settings
            'default_protein_style': self.default_protein_style.currentText(),
            'default_pigment_style': self.default_pigment_style.currentText(),
            'default_color_scheme': self.default_color_scheme.currentText(),
            'figure_width': self.figure_width_spinbox.value(),
            'figure_height': self.figure_height_spinbox.value(),
            'dpi': self.dpi_spinbox.value(),
            'font_size': self.font_size_spinbox.value(),
            'use_colorblind_friendly': self.use_colorblind_friendly.isChecked(),
            'high_contrast_mode': self.high_contrast_mode.isChecked()
        }