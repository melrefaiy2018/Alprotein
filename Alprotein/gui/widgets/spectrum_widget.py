"""
Spectrum Visualization Widget for Alprotein GUI

This widget provides visualization for absorption and fluorescence spectra with:
- Interactive plots with zoom/pan
- Peak detection and annotation
- Export capabilities
- FWHM calculation
- Experimental data loading and overlay
"""

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QCheckBox, QDoubleSpinBox, QComboBox, QTabWidget,
    QSizePolicy, QFileDialog, QMessageBox, QLineEdit, QSpinBox,
    QScrollArea, QFrame
)
from PyQt5.QtCore import Qt, pyqtSignal, QSize
from PyQt5.QtGui import QFont

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import sys


class CustomToolbar(NavigationToolbar):
    """Custom toolbar with only essential tools and bigger icons"""
    
    # Keep only essential tools: Home, Pan, Zoom, Save
    toolitems = [t for t in NavigationToolbar.toolitems if
                 t[0] in ('Home', 'Pan', 'Zoom', 'Save')]
    
    def __init__(self, canvas, parent):
        super().__init__(canvas, parent)
        # Set icon size to be clearer and more visible
        self.setIconSize(QSize(32, 32))


class SingleSpectrumPanel(QWidget):
    """
    Single spectrum panel (absorption OR fluorescence)
    """

    def __init__(self, spectrum_type="Absorption", parent=None):
        """
        Args:
            spectrum_type: "Absorption" or "Fluorescence"
        """
        super().__init__(parent)
        self.spectrum_type = spectrum_type
        self.wavelengths = None
        self.spectrum = None
        self.exp_wavelengths = None
        self.exp_spectrum = None
        self.exp_file_path = None
        
        # Plot customization settings
        self.plot_title = f"{spectrum_type} Spectrum"
        self.x_label = "Wavelength (nm)"
        self.y_label = "Normalized Intensity"
        self.x_min = None
        self.x_max = None
        self.y_min = -0.05
        self.y_max = 1.1
        self.line_width = 2.5
        self.line_style = "-"
        self.line_color = "#0000FF" if spectrum_type == "Absorption" else "#FF0000"
        self.show_fill = True
        self.show_grid = True
        
        self.setup_ui()

    def setup_ui(self):
        """Setup the UI"""
        main_layout = QHBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)
        
        # Left side: plot area
        plot_widget = QWidget()
        layout = QVBoxLayout(plot_widget)
        layout.setContentsMargins(8, 6, 8, 6)
        layout.setSpacing(6)

        # Header with title and peak info
        header_layout = QHBoxLayout()
        title = QLabel(f"{self.spectrum_type.upper()} SPECTRUM")
        title.setStyleSheet("color: #111111; font-size: 13px; font-weight: bold;")
        header_layout.addWidget(title)

        self.peak_info_label = QLabel("")
        self.peak_info_label.setStyleSheet("color: #1f7a44; font-weight: bold; font-size: 11px;")
        header_layout.addWidget(self.peak_info_label)
        header_layout.addStretch()
        layout.addLayout(header_layout)

        # Experimental data file loader - compact version
        exp_layout = QHBoxLayout()
        exp_label = QLabel("Exp Data:")
        exp_label.setStyleSheet("font-weight: bold; color: #111111; background-color: transparent; font-size: 11px;")
        exp_layout.addWidget(exp_label)

        self.exp_file_label = QLineEdit()
        self.exp_file_label.setReadOnly(True)
        self.exp_file_label.setPlaceholderText("No experimental data loaded")
        self.exp_file_label.setMaximumHeight(24)
        self.exp_file_label.setStyleSheet(
            "color: #6b6b6b; background-color: #ffffff; font-style: italic; font-size: 10px; "
            "border: 1px solid #e1e1e1; border-radius: 2px; padding: 2px 6px;"
        )
        exp_layout.addWidget(self.exp_file_label)

        self.browse_exp_btn = QPushButton("📁 Browse")
        self.browse_exp_btn.clicked.connect(self.browse_experimental_data)
        self.browse_exp_btn.setMaximumHeight(24)
        self.browse_exp_btn.setStyleSheet("""
            QPushButton {
                background-color: #111111;
                color: #ffffff;
                padding: 3px 8px;
                border-radius: 2px;
                border: 1px solid #111111;
                font-size: 10px;
            }
            QPushButton:hover {
                background-color: #000000;
            }
        """)
        exp_layout.addWidget(self.browse_exp_btn)

        self.clear_exp_btn = QPushButton("Clear")
        self.clear_exp_btn.clicked.connect(self.clear_experimental_data)
        self.clear_exp_btn.setEnabled(False)
        self.clear_exp_btn.setMaximumHeight(24)
        self.clear_exp_btn.setStyleSheet("""
            QPushButton {
                background-color: #f3f3f3;
                color: #111111;
                padding: 3px 8px;
                border-radius: 2px;
                border: 1px solid #e1e1e1;
                font-size: 10px;
            }
            QPushButton:hover {
                background-color: #ededed;
            }
            QPushButton:disabled {
                background-color: #f6f6f6;
                color: #b0b0b0;
                border: 1px solid #e9e9e9;
            }
        """)
        exp_layout.addWidget(self.clear_exp_btn)
        layout.addLayout(exp_layout)

        # Matplotlib figure with optimized size
        self.figure = Figure(figsize=(8, 5), facecolor='white', dpi=100)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.canvas.setMinimumHeight(350)

        # Toolbar with all editing capabilities - modern professional styling
        self.toolbar = CustomToolbar(self.canvas, self)
        self.toolbar.setMaximumHeight(52)
        self.toolbar.setStyleSheet("""
            QToolBar {
                background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                    stop:0 #ffffff, stop:1 #f8f9fa);
                border: none;
                border-bottom: 1px solid #e1e4e8;
                padding: 6px 10px;
                spacing: 6px;
            }
            QToolButton {
                background-color: transparent;
                border: 1px solid transparent;
                border-radius: 5px;
                padding: 8px;
                margin: 2px;
                color: #24292e;
                min-width: 36px;
                min-height: 36px;
            }
            QToolButton:hover {
                background-color: #e1e4e8;
                border: 1px solid #d1d5da;
            }
            QToolButton:pressed {
                background-color: #d1d5da;
                border: 1px solid #c1c5ca;
            }
            QToolButton:checked {
                background-color: #0366d6;
                border: 1px solid #0366d6;
                color: white;
            }
            QToolButton:disabled {
                color: #959da5;
                background-color: transparent;
            }
            QLabel {
                color: #586069;
                padding: 0px 8px;
                font-size: 11px;
            }
        """)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas, stretch=1)

        # Export buttons
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()

        self.export_csv_btn = QPushButton("📊 Export CSV")
        self.export_csv_btn.clicked.connect(self.export_csv)
        self.export_csv_btn.setStyleSheet("""
            QPushButton {
                background-color: #111111;
                color: #ffffff;
                padding: 6px 12px;
                border-radius: 2px;
                border: 1px solid #111111;
            }
            QPushButton:hover {
                background-color: #000000;
            }
        """)
        btn_layout.addWidget(self.export_csv_btn)

        self.export_plot_btn = QPushButton("📈 Export Plot")
        self.export_plot_btn.clicked.connect(self.export_plot)
        self.export_plot_btn.setStyleSheet("""
            QPushButton {
                background-color: #ffffff;
                color: #111111;
                padding: 6px 12px;
                border-radius: 2px;
                border: 1px solid #111111;
            }
            QPushButton:hover {
                background-color: #efefef;
            }
        """)
        btn_layout.addWidget(self.export_plot_btn)

        layout.addLayout(btn_layout)
        
        main_layout.addWidget(plot_widget, stretch=1)
        
        # Right side: collapsible settings panel
        self.settings_panel = self.create_settings_panel()
        main_layout.addWidget(self.settings_panel)
    
    def create_settings_panel(self):
        """Create integrated settings panel for plot customization"""
        panel = QFrame()
        panel.setMaximumWidth(280)
        panel.setStyleSheet("""
            QFrame {
                background-color: #f8f9fa;
                border-left: 1px solid #e1e4e8;
            }
        """)
        
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setStyleSheet("""
            QScrollArea { 
                border: none; 
                background-color: #f8f9fa; 
            }
            QScrollBar:vertical {
                background-color: #f8f9fa;
                width: 10px;
                border-radius: 5px;
            }
            QScrollBar::handle:vertical {
                background-color: #d1d5da;
                border-radius: 5px;
                min-height: 20px;
            }
            QScrollBar::handle:vertical:hover {
                background-color: #959da5;
            }
        """)
        
        content = QWidget()
        content.setStyleSheet("background-color: #f8f9fa;")
        panel_layout = QVBoxLayout(content)
        panel_layout.setContentsMargins(12, 12, 12, 12)
        panel_layout.setSpacing(12)
        
        # Header
        header = QLabel("⚙️ Plot Settings")
        header.setStyleSheet("""
            font-size: 13px;
            font-weight: bold;
            color: #24292e;
            padding: 4px 0px;
        """)
        panel_layout.addWidget(header)
        
        # Title section
        title_group = self.create_setting_group("Title", [
            ("Title:", QLineEdit(self.plot_title), "title_input")
        ])
        panel_layout.addWidget(title_group)
        
        # X-Axis section
        x_axis_group = self.create_setting_group("X-Axis", [
            ("Label:", QLineEdit(self.x_label), "x_label_input"),
            ("Min:", QDoubleSpinBox(), "x_min_spin"),
            ("Max:", QDoubleSpinBox(), "x_max_spin"),
        ])
        panel_layout.addWidget(x_axis_group)
        
        # Y-Axis section
        y_axis_group = self.create_setting_group("Y-Axis", [
            ("Label:", QLineEdit(self.y_label), "y_label_input"),
            ("Min:", QDoubleSpinBox(), "y_min_spin"),
            ("Max:", QDoubleSpinBox(), "y_max_spin"),
        ])
        panel_layout.addWidget(y_axis_group)
        
        # Curve section
        curve_group = self.create_setting_group("Curve Style", [
            ("Line Width:", QDoubleSpinBox(), "linewidth_spin"),
            ("Line Style:", QComboBox(), "linestyle_combo"),
            ("Color:", QLineEdit(self.line_color), "color_input"),
        ])
        panel_layout.addWidget(curve_group)
        
        # Appearance section
        appearance_group = QGroupBox("Appearance")
        appearance_group.setStyleSheet(self.get_group_style() + """
            QCheckBox {
                font-size: 11px;
                color: #111111;
                background-color: transparent;
                padding: 4px;
            }
            QCheckBox::indicator {
                width: 16px;
                height: 16px;
                border: 1px solid #e1e4e8;
                border-radius: 2px;
                background-color: #ffffff;
            }
            QCheckBox::indicator:checked {
                background-color: #0366d6;
                border: 1px solid #0366d6;
            }
        """)
        appearance_layout = QVBoxLayout(appearance_group)
        appearance_layout.setSpacing(8)
        
        self.fill_checkbox = QCheckBox("Show fill under curve")
        self.fill_checkbox.setChecked(self.show_fill)
        self.fill_checkbox.stateChanged.connect(self.on_plot_settings_changed)
        appearance_layout.addWidget(self.fill_checkbox)
        
        self.grid_checkbox = QCheckBox("Show grid")
        self.grid_checkbox.setChecked(self.show_grid)
        self.grid_checkbox.stateChanged.connect(self.on_plot_settings_changed)
        appearance_layout.addWidget(self.grid_checkbox)
        
        panel_layout.addWidget(appearance_group)
        
        # Apply button
        apply_btn = QPushButton("🔄 Apply Changes")
        apply_btn.clicked.connect(self.on_plot_settings_changed)
        apply_btn.setStyleSheet("""
            QPushButton {
                background-color: #111111;
                color: #ffffff;
                padding: 8px 16px;
                border-radius: 2px;
                border: 1px solid #111111;
                font-weight: bold;
                font-size: 11px;
            }
            QPushButton:hover {
                background-color: #000000;
            }
            QPushButton:pressed {
                background-color: #000000;
            }
        """)
        panel_layout.addWidget(apply_btn)
        
        panel_layout.addStretch()
        
        scroll.setWidget(content)
        
        panel_vlayout = QVBoxLayout(panel)
        panel_vlayout.setContentsMargins(0, 0, 0, 0)
        panel_vlayout.addWidget(scroll)
        
        # Configure spinboxes
        self.x_min_spin.setRange(-10000, 10000)
        self.x_min_spin.setValue(self.x_min if self.x_min is not None else 0)
        self.x_min_spin.setSpecialValueText("Auto")
        self.x_min_spin.valueChanged.connect(lambda v: setattr(self, 'x_min', v if v != 0 else None))
        
        self.x_max_spin.setRange(-10000, 10000)
        self.x_max_spin.setValue(self.x_max if self.x_max is not None else 0)
        self.x_max_spin.setSpecialValueText("Auto")
        self.x_max_spin.valueChanged.connect(lambda v: setattr(self, 'x_max', v if v != 0 else None))
        
        self.y_min_spin.setRange(-100, 100)
        self.y_min_spin.setValue(self.y_min)
        self.y_min_spin.setSingleStep(0.1)
        self.y_min_spin.setDecimals(2)
        self.y_min_spin.valueChanged.connect(lambda v: setattr(self, 'y_min', v))
        
        self.y_max_spin.setRange(-100, 100)
        self.y_max_spin.setValue(self.y_max)
        self.y_max_spin.setSingleStep(0.1)
        self.y_max_spin.setDecimals(2)
        self.y_max_spin.valueChanged.connect(lambda v: setattr(self, 'y_max', v))
        
        self.linewidth_spin.setRange(0.5, 10)
        self.linewidth_spin.setValue(self.line_width)
        self.linewidth_spin.setSingleStep(0.5)
        self.linewidth_spin.setDecimals(1)
        self.linewidth_spin.valueChanged.connect(lambda v: setattr(self, 'line_width', v))
        
        self.linestyle_combo.addItems(["-", "--", "-.", ":"])
        self.linestyle_combo.setCurrentText(self.line_style)
        self.linestyle_combo.currentTextChanged.connect(lambda v: setattr(self, 'line_style', v))
        
        self.title_input.textChanged.connect(lambda v: setattr(self, 'plot_title', v))
        self.x_label_input.textChanged.connect(lambda v: setattr(self, 'x_label', v))
        self.y_label_input.textChanged.connect(lambda v: setattr(self, 'y_label', v))
        self.color_input.textChanged.connect(lambda v: setattr(self, 'line_color', v))
        
        return panel
    
    def create_setting_group(self, title, widgets):
        """Create a settings group box"""
        group = QGroupBox(title)
        group.setStyleSheet(self.get_group_style())
        layout = QVBoxLayout(group)
        layout.setSpacing(8)
        
        for label_text, widget, attr_name in widgets:
            row = QHBoxLayout()
            label = QLabel(label_text)
            label.setStyleSheet("""
                color: #111111; 
                font-size: 11px;
                font-weight: normal;
                background-color: transparent;
            """)
            label.setMinimumWidth(70)
            row.addWidget(label)
            
            widget.setStyleSheet(self.get_input_style())
            setattr(self, attr_name, widget)
            row.addWidget(widget, 1)
            layout.addLayout(row)
        
        return group
    
    def get_group_style(self):
        """Get group box stylesheet"""
        return """
            QGroupBox {
                font-weight: bold;
                font-size: 11px;
                color: #111111;
                border: 1px solid #e1e4e8;
                border-radius: 3px;
                margin-top: 8px;
                padding-top: 8px;
                background-color: #ffffff;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                subcontrol-position: top left;
                left: 8px;
                padding: 0 4px 0 4px;
                color: #111111;
                background-color: #ffffff;
            }
        """
    
    def get_input_style(self):
        """Get input widget stylesheet"""
        return """
            QLineEdit, QDoubleSpinBox, QSpinBox, QComboBox {
                background-color: #ffffff;
                border: 1px solid #e1e4e8;
                border-radius: 2px;
                padding: 4px 6px;
                font-size: 11px;
                color: #111111;
            }
            QLineEdit:focus, QDoubleSpinBox:focus, QSpinBox:focus, QComboBox:focus {
                border: 1px solid #111111;
                background-color: #ffffff;
            }
            QDoubleSpinBox::up-button, QDoubleSpinBox::down-button,
            QSpinBox::up-button, QSpinBox::down-button {
                background-color: #f3f3f3;
                border: 1px solid #e1e4e8;
            }
            QDoubleSpinBox::up-button:hover, QDoubleSpinBox::down-button:hover,
            QSpinBox::up-button:hover, QSpinBox::down-button:hover {
                background-color: #e1e4e8;
            }
            QComboBox::drop-down {
                border: none;
                background-color: #f3f3f3;
            }
            QComboBox::down-arrow {
                image: none;
                border-left: 4px solid transparent;
                border-right: 4px solid transparent;
                border-top: 5px solid #111111;
                width: 0;
                height: 0;
            }
            QComboBox QAbstractItemView {
                background-color: #ffffff;
                color: #111111;
                border: 1px solid #e1e4e8;
                selection-background-color: #f3f3f3;
                selection-color: #111111;
            }
            QCheckBox {
                font-size: 11px;
                color: #111111;
                background-color: transparent;
            }
            QCheckBox::indicator {
                width: 16px;
                height: 16px;
                border: 1px solid #e1e4e8;
                border-radius: 2px;
                background-color: #ffffff;
            }
            QCheckBox::indicator:checked {
                background-color: #0366d6;
                border: 1px solid #0366d6;
            }
        """
    
    def on_plot_settings_changed(self):
        """Handle plot settings changes"""
        self.show_fill = self.fill_checkbox.isChecked()
        self.show_grid = self.grid_checkbox.isChecked()
        self.update_plot()

    def browse_experimental_data(self):
        """Browse for experimental data file"""
        options = QFileDialog.DontUseNativeDialog if sys.platform == 'darwin' else QFileDialog.Options()

        file_path, _ = QFileDialog.getOpenFileName(
            self,
            f"Load Experimental {self.spectrum_type} Data",
            "",
            "NumPy Files (*.npy);;CSV Files (*.csv);;All Files (*)",
            options=options
        )

        if file_path:
            self.load_experimental_data(file_path)

    def load_experimental_data(self, file_path):
        """Load experimental data from file"""
        try:
            if file_path.endswith('.npy'):
                data = np.load(file_path)
                if data.ndim == 2 and data.shape[1] >= 2:
                    self.exp_wavelengths = data[:, 0]
                    self.exp_spectrum = data[:, 1]
                    # Normalize
                    if np.max(self.exp_spectrum) > 0:
                        self.exp_spectrum = self.exp_spectrum / np.max(self.exp_spectrum)
                else:
                    raise ValueError(f"Unexpected data shape: {data.shape}")
            elif file_path.endswith('.csv'):
                data = np.loadtxt(file_path, delimiter=',', skiprows=1)
                if data.ndim == 2 and data.shape[1] >= 2:
                    self.exp_wavelengths = data[:, 0]
                    self.exp_spectrum = data[:, 1]
                    # Normalize
                    if np.max(self.exp_spectrum) > 0:
                        self.exp_spectrum = self.exp_spectrum / np.max(self.exp_spectrum)
                else:
                    raise ValueError(f"Unexpected data shape: {data.shape}")
            else:
                raise ValueError("Unsupported file format")

            self.exp_file_path = file_path
            self.exp_file_label.setText(Path(file_path).name)
            self.exp_file_label.setStyleSheet(
                "color: #111111; background-color: #ffffff; "
                "border: 1px solid #e1e1e1; border-radius: 2px; padding: 2px 6px;"
            )
            self.clear_exp_btn.setEnabled(True)

            # Update plot
            self.update_plot()

            QMessageBox.information(self, "Success",
                                  f"Loaded experimental {self.spectrum_type.lower()} data:\n{Path(file_path).name}")
        except Exception as e:
            QMessageBox.critical(self, "Error",
                               f"Failed to load experimental data:\n{str(e)}")

    def clear_experimental_data(self):
        """Clear experimental data"""
        self.exp_wavelengths = None
        self.exp_spectrum = None
        self.exp_file_path = None
        self.exp_file_label.setText("")
        self.exp_file_label.setPlaceholderText("No experimental data loaded")
        self.exp_file_label.setStyleSheet(
            "color: #6b6b6b; background-color: #ffffff; font-style: italic; "
            "border: 1px solid #e1e1e1; border-radius: 2px; padding: 2px 6px;"
        )
        self.clear_exp_btn.setEnabled(False)
        self.update_plot()

    def update_spectrum(self, wavelengths: np.ndarray, spectrum: np.ndarray):
        """Update with simulated spectrum data"""
        print(f"[DEBUG] {self.spectrum_type} panel update_spectrum called")
        print(f"  wavelengths: {wavelengths.shape if wavelengths is not None else None}")
        print(f"  spectrum: {spectrum.shape if spectrum is not None else None}")
        self.wavelengths = wavelengths
        self.spectrum = spectrum
        self.update_plot()

    def update_plot(self):
        """Update the plot"""
        print(f"[DEBUG] {self.spectrum_type} panel update_plot called")
        if self.wavelengths is None or self.spectrum is None:
            print(f"  Skipping plot update - wavelengths or spectrum is None")
            return
        print(f"  Proceeding with plot update...")

        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.set_facecolor('white')

        # Plot experimental data if available
        if self.exp_wavelengths is not None and self.exp_spectrum is not None:
            ax.plot(self.exp_wavelengths, self.exp_spectrum, 'k-',
                   linewidth=2, alpha=0.7, label='Experimental')

        # Plot simulated spectrum with customization
        ax.plot(self.wavelengths, self.spectrum, 
               linestyle=self.line_style,
               color=self.line_color,
               linewidth=self.line_width, 
               label=f'Simulated (Renger)')
        
        if self.show_fill:
            ax.fill_between(self.wavelengths, self.spectrum, alpha=0.2, color=self.line_color)

        # Find peak (for info display, but don't plot the star)
        peak_idx = np.argmax(self.spectrum)
        peak_wavelength = self.wavelengths[peak_idx]
        peak_intensity = self.spectrum[peak_idx]

        # Calculate FWHM
        fwhm = self.calculate_fwhm(self.wavelengths, self.spectrum)

        # Update peak info label
        self.peak_info_label.setText(
            f"Peak: {peak_wavelength:.1f} nm ({1e7/peak_wavelength:.1f} cm⁻¹) | FWHM: {fwhm:.1f} nm"
        )

        # Labels and formatting with customization
        ax.set_xlabel(self.x_label, fontsize=12, color='black', fontweight='bold')
        ax.set_ylabel(self.y_label, fontsize=12, color='black', fontweight='bold')
        ax.set_title(self.plot_title,
                    fontsize=13, color='#111111', fontweight='bold', pad=10)
        ax.legend(loc='upper right', framealpha=0.9, fontsize=10)
        
        if self.show_grid:
            ax.grid(True, alpha=0.25, linestyle='--', linewidth=0.5)

        # Set axis limits with customization
        if len(self.wavelengths) > 0:
            x_min = self.x_min if self.x_min is not None else self.wavelengths.min()
            x_max = self.x_max if self.x_max is not None else self.wavelengths.max()
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(self.y_min, self.y_max)

        # Set tick colors
        ax.tick_params(colors='black', labelsize=10)
        for spine in ax.spines.values():
            spine.set_edgecolor('#e1e1e1')
            spine.set_linewidth(1.5)

        self.figure.tight_layout()
        self.canvas.draw()

    def calculate_fwhm(self, wavelengths: np.ndarray, spectrum: np.ndarray) -> float:
        """Calculate Full Width at Half Maximum"""
        peak_idx = np.argmax(spectrum)
        peak_value = spectrum[peak_idx]
        half_max = peak_value / 2

        # Find left and right points at half maximum
        left_idx = np.where(spectrum[:peak_idx] <= half_max)[0]
        right_idx = np.where(spectrum[peak_idx:] <= half_max)[0]

        if len(left_idx) > 0 and len(right_idx) > 0:
            left_wl = wavelengths[left_idx[-1]]
            right_wl = wavelengths[peak_idx + right_idx[0]]
            fwhm = right_wl - left_wl
        else:
            fwhm = 0.0

        return fwhm

    def export_csv(self):
        """Export spectrum data to CSV"""
        if self.wavelengths is None or self.spectrum is None:
            QMessageBox.warning(self, "No Data", "No spectrum data to export")
            return

        options = QFileDialog.DontUseNativeDialog if sys.platform == 'darwin' else QFileDialog.Options()

        default_name = f"{self.spectrum_type.lower()}_spectrum.csv"
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            f"Export {self.spectrum_type} Spectrum",
            default_name,
            "CSV Files (*.csv);;All Files (*)",
            options=options
        )

        if file_path:
            try:
                import csv
                with open(file_path, 'w', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow(['wavelength_nm', self.spectrum_type.lower()])
                    for wl, intensity in zip(self.wavelengths, self.spectrum):
                        writer.writerow([f"{wl:.2f}", f"{intensity:.6f}"])

                QMessageBox.information(self, "Success",
                                      f"{self.spectrum_type} spectrum exported to:\n{file_path}")
            except Exception as e:
                QMessageBox.critical(self, "Error",
                                   f"Failed to export:\n{str(e)}")

    def export_plot(self):
        """Export plot as image"""
        options = QFileDialog.DontUseNativeDialog if sys.platform == 'darwin' else QFileDialog.Options()

        default_name = f"{self.spectrum_type.lower()}_spectrum.png"
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            f"Export {self.spectrum_type} Plot",
            default_name,
            "PNG Files (*.png);;PDF Files (*.pdf);;All Files (*)",
            options=options
        )

        if file_path:
            self.figure.savefig(file_path, dpi=300, bbox_inches='tight',
                              facecolor='white', edgecolor='none')
            QMessageBox.information(self, "Success",
                                  f"Plot exported to:\n{file_path}")

    def clear(self):
        """Clear all data"""
        self.wavelengths = None
        self.spectrum = None
        self.clear_experimental_data()
        self.figure.clear()
        self.canvas.draw()
        self.peak_info_label.setText("")


class SpectrumPlotWidget(QWidget):
    """
    Widget for plotting both absorption and fluorescence spectra
    Uses tab widget to separate absorption and fluorescence panels
    """

    export_requested = pyqtSignal()
    run_spectra_requested = pyqtSignal()
    parameters_changed = pyqtSignal(dict)

    def __init__(self, parent=None):
        super().__init__(parent)
        # Store data for both spectra
        self.wavelengths_abs = None
        self.absorption = None
        self.wavelengths_fl = None
        self.fluorescence = None
        # Legacy data for backward compatibility
        self.spectrum_data = None
        self.wavelengths = None
        self.intensities = None
        self.stick_wavelengths = None
        self.stick_intensities = None
        self.setup_ui()

    def setup_ui(self):
        """Setup the UI with tabs for absorption and fluorescence"""
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(8, 6, 8, 6)
        main_layout.setSpacing(6)

        header_layout = QHBoxLayout()
        header_label = QLabel("SPECTRA")
        header_label.setStyleSheet("color: #111111; font-size: 13px; font-weight: bold;")
        header_layout.addWidget(header_label)
        header_layout.addStretch()

        self.run_spectra_btn = QPushButton("▶ Run Spectra")
        self.run_spectra_btn.clicked.connect(self.run_spectra_requested.emit)
        self.run_spectra_btn.setEnabled(False)
        self.run_spectra_btn.setStyleSheet("""
            QPushButton {
                background-color: #111111;
                color: #ffffff;
                border: 1px solid #111111;
                padding: 5px 12px;
                border-radius: 2px;
                font-size: 11px;
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
        header_layout.addWidget(self.run_spectra_btn)
        main_layout.addLayout(header_layout)

        # Spectra settings - compact horizontal layout
        settings_group = QGroupBox("SPECTRA SETTINGS")
        settings_group.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                border: 1px solid #e1e1e1;
                border-radius: 2px;
                margin-top: 4px;
                padding-top: 6px;
                background-color: #ffffff;
                color: #111111;
                font-size: 11px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 8px;
                padding: 0 4px;
            }
            QLabel {
                color: #111111;
                background-color: transparent;
                font-weight: normal;
                font-size: 11px;
            }
            QSpinBox {
                color: #111111;
                background-color: #ffffff;
                border: 1px solid #e1e1e1;
                border-radius: 2px;
                padding: 3px;
                font-size: 11px;
                max-height: 24px;
            }
            QSpinBox::up-button, QSpinBox::down-button {
                background-color: #f3f3f3;
                border: 1px solid #e1e1e1;
            }
            QSpinBox::up-button:hover, QSpinBox::down-button:hover {
                background-color: #e1e1e1;
            }
            QCheckBox {
                color: #111111;
                background-color: transparent;
                font-size: 11px;
            }
        """)
        settings_layout = QHBoxLayout(settings_group)
        settings_layout.setSpacing(8)
        settings_layout.setContentsMargins(8, 6, 8, 6)

        # Temperature
        temp_label = QLabel("Temperature")
        settings_layout.addWidget(temp_label)
        self.temperature = QSpinBox()
        self.temperature.setRange(0, 500)
        self.temperature.setValue(300)
        self.temperature.setSuffix(" K")
        self.temperature.setMaximumWidth(100)
        self.temperature.valueChanged.connect(self.on_parameters_changed)
        settings_layout.addWidget(self.temperature)

        # N Ensemble
        ensemble_label = QLabel("N Ensemble")
        settings_layout.addWidget(ensemble_label)
        self.n_ensemble = QSpinBox()
        self.n_ensemble.setRange(10, 500)
        self.n_ensemble.setValue(100)
        self.n_ensemble.setMaximumWidth(80)
        self.n_ensemble.valueChanged.connect(self.on_parameters_changed)
        settings_layout.addWidget(self.n_ensemble)

        # Checkboxes
        self.include_vibronic = QCheckBox("Include vibronic (0-1) transitions")
        self.include_vibronic.setChecked(True)
        self.include_vibronic.stateChanged.connect(self.on_parameters_changed)
        settings_layout.addWidget(self.include_vibronic)

        self.use_ensemble = QCheckBox("Use ensemble averaging")
        self.use_ensemble.setChecked(True)
        self.use_ensemble.stateChanged.connect(self.on_parameters_changed)
        settings_layout.addWidget(self.use_ensemble)

        settings_layout.addStretch()

        main_layout.addWidget(settings_group)

        # Create tab widget for absorption and fluorescence
        self.tab_widget = QTabWidget()
        self.tab_widget.setStyleSheet("""
            QTabWidget::pane {
                background-color: #ffffff;
                border: 1px solid #e1e1e1;
                border-radius: 2px;
            }
            QTabBar {
                qproperty-expanding: 0;
            }
            QTabBar::tab {
                color: #111111;
                background-color: #f3f3f3;
                padding: 6px 20px;
                min-width: 90px;
                margin: 2px;
                border: 1px solid #e1e1e1;
                border-bottom: none;
                border-top-left-radius: 2px;
                border-top-right-radius: 2px;
                font-size: 12px;
            }
            QTabBar::tab:selected {
                background-color: #ffffff;
                color: #111111;
                font-weight: bold;
                border-bottom: 2px solid #111111;
            }
            QTabBar::tab:hover {
                background-color: #ededed;
                color: #111111;
            }
        """)

        # Create panels for absorption and fluorescence
        self.absorption_panel = SingleSpectrumPanel("Absorption")
        self.fluorescence_panel = SingleSpectrumPanel("Fluorescence")

        # Add panels to tab widget
        self.tab_widget.addTab(self.absorption_panel, "📊 Absorption")
        self.tab_widget.addTab(self.fluorescence_panel, "🌟 Fluorescence")

        main_layout.addWidget(self.tab_widget, stretch=1)

    def get_parameters(self) -> Dict[str, Any]:
        """Return spectra-specific parameters"""
        return {
            'temperature': self.temperature.value(),
            'n_ensemble': self.n_ensemble.value(),
            'include_vibronic': self.include_vibronic.isChecked(),
            'use_ensemble': self.use_ensemble.isChecked(),
        }

    def on_parameters_changed(self):
        """Emit parameter changes for spectra settings"""
        self.parameters_changed.emit(self.get_parameters())

    def set_run_enabled(self, enabled: bool):
        """Enable or disable the Run Spectra button"""
        self.run_spectra_btn.setEnabled(enabled)

    def update_spectra(self, wavelengths_abs: np.ndarray, absorption: np.ndarray,
                      wavelengths_fl: np.ndarray, fluorescence: np.ndarray):
        """
        Update both absorption and fluorescence spectra

        Args:
            wavelengths_abs: Wavelength axis for absorption (nm)
            absorption: Absorption spectrum (normalized)
            wavelengths_fl: Wavelength axis for fluorescence (nm)
            fluorescence: Fluorescence spectrum (normalized)
        """
        # Debug logging
        print(f"[DEBUG] update_spectra called")
        print(f"  Absorption: wavelengths shape={wavelengths_abs.shape if wavelengths_abs is not None else None}, "
              f"spectrum shape={absorption.shape if absorption is not None else None}")
        print(f"  Fluorescence: wavelengths shape={wavelengths_fl.shape if wavelengths_fl is not None else None}, "
              f"spectrum shape={fluorescence.shape if fluorescence is not None else None}")

        self.wavelengths_abs = wavelengths_abs
        self.absorption = absorption
        self.wavelengths_fl = wavelengths_fl
        self.fluorescence = fluorescence

        # Update individual panels
        print(f"  Updating absorption panel...")
        self.absorption_panel.update_spectrum(wavelengths_abs, absorption)
        print(f"  Updating fluorescence panel...")
        self.fluorescence_panel.update_spectrum(wavelengths_fl, fluorescence)
        print(f"  Both panels updated")

    def update_spectrum_data(self, eigenvalues: np.ndarray, oscillator_strengths: np.ndarray):
        """
        Legacy method for backward compatibility
        Updates only absorption spectrum from eigenvalues

        Args:
            eigenvalues: Array of energy eigenvalues in cm^-1
            oscillator_strengths: Array of oscillator strengths
        """
        # Store for legacy compatibility
        self.stick_wavelengths = 1e7 / eigenvalues
        self.stick_intensities = oscillator_strengths

        # Generate a simple broadened spectrum for display
        wavelengths = np.linspace(600, 750, 1000)
        spectrum = np.zeros_like(wavelengths)

        # Simple Gaussian broadening
        for wl, intensity in zip(self.stick_wavelengths, self.stick_intensities):
            sigma = 10.0  # nm
            spectrum += intensity * np.exp(-((wavelengths - wl) ** 2) / (2 * sigma ** 2))

        # Normalize
        if spectrum.max() > 0:
            spectrum /= spectrum.max()

        # Update only absorption panel
        self.absorption_panel.update_spectrum(wavelengths, spectrum)

    def update_spectrum_with_components(self, wavelengths: np.ndarray, spectrum: np.ndarray,
                                       A_00: Optional[np.ndarray] = None,
                                       A_01: Optional[np.ndarray] = None):
        """
        Legacy method for backward compatibility with Renger lineshape
        Updates only absorption spectrum

        Args:
            wavelengths: Wavelength axis (nm)
            spectrum: Total absorption spectrum (normalized)
            A_00: 0-0 (zero-phonon) component
            A_01: 0-1 (vibronic) component
        """
        self.wavelengths = wavelengths
        self.spectrum = spectrum
        self.A_00 = A_00
        self.A_01 = A_01

        # Update only absorption panel
        self.absorption_panel.update_spectrum(wavelengths, spectrum)

    def load_experimental_absorption(self, file_path: str):
        """Load experimental absorption data"""
        self.absorption_panel.load_experimental_data(file_path)

    def load_experimental_fluorescence(self, file_path: str):
        """Load experimental fluorescence data"""
        self.fluorescence_panel.load_experimental_data(file_path)

    def clear(self):
        """Clear both panels"""
        self.absorption_panel.clear()
        self.fluorescence_panel.clear()
        self.wavelengths_abs = None
        self.absorption = None
        self.wavelengths_fl = None
        self.fluorescence = None
        # Clear legacy data
        self.spectrum_data = None
        self.wavelengths = None
        self.intensities = None
        self.stick_wavelengths = None
        self.stick_intensities = None
