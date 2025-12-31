"""
Control panel for 3D visualization settings
"""

from typing import Dict, Any
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton, 
                            QLabel, QComboBox, QSlider, QCheckBox, QGroupBox,
                            QListWidget, QListWidgetItem)
from PyQt5.QtCore import Qt, pyqtSignal


class ControlPanel3D(QWidget):
    """Control panel for 3D visualization settings"""
    
    view_changed = pyqtSignal(dict)  # settings
    pigment_selected = pyqtSignal(str)  # pigment_id
    
    def __init__(self):
        super().__init__()
        self.setup_ui()
        self.connect_signals()
    
    def setup_ui(self):
        """Setup control panel UI"""
        layout = QVBoxLayout(self)
        
        # Visualization Style Group
        style_group = QGroupBox("🎨 Visualization Style")
        style_layout = QVBoxLayout(style_group)
        
        # Protein representation
        style_layout.addWidget(QLabel("Protein Representation:"))
        self.protein_style = QComboBox()
        self.protein_style.addItems([
            "Cartoon", "Surface", "Line", "Stick", "Sphere"
        ])
        self.protein_style.setCurrentText("Cartoon")
        style_layout.addWidget(self.protein_style)
        
        # Pigment representation
        style_layout.addWidget(QLabel("Pigment Representation:"))
        self.pigment_style = QComboBox()
        self.pigment_style.addItems([
            "Stick", "Ball & Stick", "Sphere", "Line"
        ])
        self.pigment_style.setCurrentText("Stick")
        style_layout.addWidget(self.pigment_style)
        
        # Color scheme
        style_layout.addWidget(QLabel("Color Scheme:"))
        self.color_scheme = QComboBox()
        self.color_scheme.addItems([
            "By Element", "By Residue", "By Chain", "By Energy", "By Contribution"
        ])
        self.color_scheme.setCurrentText("By Element")
        style_layout.addWidget(self.color_scheme)
        
        layout.addWidget(style_group)
        
        # Display Options Group
        display_group = QGroupBox("👁️ Display Options")
        display_layout = QVBoxLayout(display_group)
        
        self.show_protein = QCheckBox("Show Protein")
        self.show_protein.setChecked(True)
        display_layout.addWidget(self.show_protein)
        
        self.show_pigments = QCheckBox("Show Pigments")
        self.show_pigments.setChecked(True)
        display_layout.addWidget(self.show_pigments)
        
        self.show_water = QCheckBox("Show Water")
        self.show_water.setChecked(False)
        display_layout.addWidget(self.show_water)
        
        self.show_ions = QCheckBox("Show Ions")
        self.show_ions.setChecked(True)
        display_layout.addWidget(self.show_ions)
        
        layout.addWidget(display_group)
        
        # Pigment Selection Group
        pigment_group = QGroupBox("🧬 Pigments")
        pigment_layout = QVBoxLayout(pigment_group)
        
        self.pigment_list = QListWidget()
        self.pigment_list.setMaximumHeight(150)
        pigment_layout.addWidget(self.pigment_list)
        
        # Pigment control buttons
        button_layout = QHBoxLayout()
        self.focus_button = QPushButton("Focus")
        self.focus_button.setEnabled(False)
        self.hide_button = QPushButton("Hide")
        self.hide_button.setEnabled(False)
        button_layout.addWidget(self.focus_button)
        button_layout.addWidget(self.hide_button)
        pigment_layout.addLayout(button_layout)
        
        layout.addWidget(pigment_group)
        
        layout.addStretch()
    
    def connect_signals(self):
        """Connect control signals"""
        # Style changes
        self.protein_style.currentTextChanged.connect(self.emit_view_changed)
        self.pigment_style.currentTextChanged.connect(self.emit_view_changed)
        self.color_scheme.currentTextChanged.connect(self.emit_view_changed)
        
        # Display options
        self.show_protein.toggled.connect(self.emit_view_changed)
        self.show_pigments.toggled.connect(self.emit_view_changed)
        self.show_water.toggled.connect(self.emit_view_changed)
        self.show_ions.toggled.connect(self.emit_view_changed)
        
        # Pigment list
        self.pigment_list.itemSelectionChanged.connect(self.on_pigment_selection_changed)
        self.focus_button.clicked.connect(self.on_focus_pigment)
        self.hide_button.clicked.connect(self.on_hide_pigment)
    
    def emit_view_changed(self):
        """Emit view changed signal with current settings"""
        settings = {
            'protein_style': self.protein_style.currentText(),
            'pigment_style': self.pigment_style.currentText(),
            'color_scheme': self.color_scheme.currentText(),
            'show_protein': self.show_protein.isChecked(),
            'show_pigments': self.show_pigments.isChecked(),
            'show_water': self.show_water.isChecked(),
            'show_ions': self.show_ions.isChecked()
        }
        self.view_changed.emit(settings)
    
    def on_pigment_selection_changed(self):
        """Handle pigment selection change"""
        items = self.pigment_list.selectedItems()
        has_selection = len(items) > 0
        self.focus_button.setEnabled(has_selection)
        self.hide_button.setEnabled(has_selection)
        
        if has_selection:
            pigment_id = items[0].text().split(':')[0]  # Extract ID from "ID: Name" format
            self.pigment_selected.emit(pigment_id)
    
    def on_focus_pigment(self):
        """Focus on selected pigment"""
        items = self.pigment_list.selectedItems()
        if items:
            pigment_id = items[0].text().split(':')[0]
            self.view_changed.emit({'action': 'focus_pigment', 'pigment_id': pigment_id})
    
    def on_hide_pigment(self):
        """Hide selected pigment"""
        items = self.pigment_list.selectedItems()
        if items:
            pigment_id = items[0].text().split(':')[0]
            self.view_changed.emit({'action': 'hide_pigment', 'pigment_id': pigment_id})
    
    def update_pigment_list(self, pigments: Dict[str, Any]):
        """Update the pigment list"""
        self.pigment_list.clear()
        for pigment_id, pigment in pigments.items():
            item_text = f"{pigment_id}: {pigment.get_resname()}"
            item = QListWidgetItem(item_text)
            self.pigment_list.addItem(item)
