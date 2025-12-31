"""
Export Dialog for results export options
"""

import os
from typing import Dict, Any
from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QGroupBox,
                            QCheckBox, QComboBox, QLineEdit, QPushButton,
                            QFileDialog, QLabel, QSpinBox, QTextEdit, 
                            QDialogButtonBox, QFormLayout, QFrame)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont


class ExportDialog(QDialog):
    """
    Dialog for configuring export options
    """
    
    def __init__(self, results: Dict[str, Any], parent=None):
        super().__init__(parent)
        self.results = results
        self.export_settings = {}
        
        self.setWindowTitle("Export Results")
        self.setModal(True)
        self.setMinimumSize(500, 600)
        
        self.setup_ui()
        self.load_default_settings()
    
    def setup_ui(self):
        """Setup the dialog UI"""
        layout = QVBoxLayout(self)
        
        # Title
        title_label = QLabel("📁 Export Calculation Results")
        title_label.setStyleSheet("font-size: 16px; font-weight: bold; margin-bottom: 10px;")
        layout.addWidget(title_label)
        
        # Output directory selection
        self.create_output_group(layout)
        
        # Export content options
        self.create_content_group(layout)
        
        # Format options
        self.create_format_group(layout)
        
        # Advanced options
        self.create_advanced_group(layout)
        
        # Preview
        self.create_preview_group(layout)
        
        # Buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)
    
    def create_output_group(self, parent_layout):
        """Create output directory selection group"""
        group = QGroupBox("📂 Output Location")
        layout = QFormLayout(group)
        
        # Directory selection
        dir_layout = QHBoxLayout()
        self.output_dir_edit = QLineEdit()
        self.output_dir_edit.setPlaceholderText("Select output directory...")
        dir_layout.addWidget(self.output_dir_edit)
        
        browse_button = QPushButton("Browse...")
        browse_button.clicked.connect(self.browse_output_dir)
        dir_layout.addWidget(browse_button)
        
        layout.addRow("Directory:", dir_layout)
        
        # Base filename
        self.base_filename_edit = QLineEdit("alprotein_results")
        layout.addRow("Base Filename:", self.base_filename_edit)
        
        parent_layout.addWidget(group)
    
    def create_content_group(self, parent_layout):
        """Create export content selection group"""
        group = QGroupBox("📋 Content to Export")
        layout = QVBoxLayout(group)
        
        # Available content based on results
        self.export_site_energies = QCheckBox("Site Energies")
        self.export_site_energies.setChecked(True)
        self.export_site_energies.setEnabled('site_energies' in self.results)
        layout.addWidget(self.export_site_energies)
        
        self.export_cdc = QCheckBox("CDC Analysis")
        self.export_cdc.setChecked(True)
        self.export_cdc.setEnabled('cdc_analysis' in self.results)
        layout.addWidget(self.export_cdc)
        
        self.export_visualization = QCheckBox("3D Visualization Files")
        self.export_visualization.setChecked(True)
        layout.addWidget(self.export_visualization)
        
        self.export_summary = QCheckBox("Summary Report")
        self.export_summary.setChecked(True)
        layout.addWidget(self.export_summary)
        
        parent_layout.addWidget(group)
    
    def create_format_group(self, parent_layout):
        """Create format options group"""
        group = QGroupBox("💾 Export Formats")
        layout = QVBoxLayout(group)
        
        # File formats
        self.format_csv = QCheckBox("CSV (Comma-separated values)")
        self.format_csv.setChecked(True)
        layout.addWidget(self.format_csv)
        
        self.format_json = QCheckBox("JSON (JavaScript Object Notation)")
        self.format_json.setChecked(True)
        layout.addWidget(self.format_json)
        
        self.format_txt = QCheckBox("TXT (Plain text reports)")
        self.format_txt.setChecked(True)
        layout.addWidget(self.format_txt)
        
        # Visualization formats
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        layout.addWidget(separator)
        
        viz_label = QLabel("Visualization Formats:")
        viz_label.setStyleSheet("font-weight: bold;")
        layout.addWidget(viz_label)
        
        self.format_png = QCheckBox("PNG (High-resolution images)")
        self.format_png.setChecked(True)
        layout.addWidget(self.format_png)
        
        self.format_html = QCheckBox("HTML (Interactive 3D viewer)")
        self.format_html.setChecked(True)
        layout.addWidget(self.format_html)
        
        self.format_pymol = QCheckBox("PyMOL Script (.pml)")
        self.format_pymol.setChecked(True)
        layout.addWidget(self.format_pymol)
        
        parent_layout.addWidget(group)
    
    def create_advanced_group(self, parent_layout):
        """Create advanced options group"""
        group = QGroupBox("⚙️ Advanced Options")
        layout = QFormLayout(group)
        
        # CDC cutoff
        self.cdc_cutoff_spin = QSpinBox()
        self.cdc_cutoff_spin.setRange(1, 1000)
        self.cdc_cutoff_spin.setValue(20)
        self.cdc_cutoff_spin.setSuffix(" cm⁻¹")
        layout.addRow("CDC Cutoff:", self.cdc_cutoff_spin)
        
        # Decimal places
        self.decimal_places_spin = QSpinBox()
        self.decimal_places_spin.setRange(0, 6)
        self.decimal_places_spin.setValue(1)
        layout.addRow("Decimal Places:", self.decimal_places_spin)
        
        # Include metadata
        self.include_metadata = QCheckBox("Include calculation metadata")
        self.include_metadata.setChecked(True)
        layout.addRow("", self.include_metadata)
        
        # Compress output
        self.compress_output = QCheckBox("Create ZIP archive")
        self.compress_output.setChecked(False)
        layout.addRow("", self.compress_output)
        
        parent_layout.addWidget(group)
    
    def create_preview_group(self, parent_layout):
        """Create export preview group"""
        group = QGroupBox("👁️ Export Preview")
        layout = QVBoxLayout(group)
        
        self.preview_text = QTextEdit()
        self.preview_text.setMaximumHeight(150)
        self.preview_text.setFont(QFont("Courier New", 9))
        self.preview_text.setReadOnly(True)
        layout.addWidget(self.preview_text)
        
        # Update preview when settings change
        self.output_dir_edit.textChanged.connect(self.update_preview)
        self.base_filename_edit.textChanged.connect(self.update_preview)
        self.export_site_energies.toggled.connect(self.update_preview)
        self.export_cdc.toggled.connect(self.update_preview)
        self.format_csv.toggled.connect(self.update_preview)
        self.format_json.toggled.connect(self.update_preview)
        
        parent_layout.addWidget(group)
    
    def browse_output_dir(self):
        """Open directory browser"""
        directory = QFileDialog.getExistingDirectory(
            self, "Select Output Directory"
        )
        if directory:
            self.output_dir_edit.setText(directory)
    
    def load_default_settings(self):
        """Load default export settings"""
        # Set default output directory to user's Documents folder
        default_dir = os.path.expanduser("~/Documents/Alprotein_Export")
        self.output_dir_edit.setText(default_dir)
        
        self.update_preview()
    
    def update_preview(self):
        """Update the export preview"""
        output_dir = self.output_dir_edit.text()
        base_filename = self.base_filename_edit.text()
        
        if not output_dir or not base_filename:
            self.preview_text.setText("Please specify output directory and filename.")
            return
        
        preview_lines = [f"Export Directory: {output_dir}", ""]
        
        # Preview files that will be created
        if self.export_site_energies.isChecked() and self.export_site_energies.isEnabled():
            preview_lines.append("Site Energies:")
            if self.format_csv.isChecked():
                preview_lines.append(f"  • {base_filename}_site_energies.csv")
            if self.format_json.isChecked():
                preview_lines.append(f"  • {base_filename}_site_energies.json")
            preview_lines.append("")
        
        if self.export_cdc.isChecked() and self.export_cdc.isEnabled():
            preview_lines.append("CDC Analysis:")
            if self.format_json.isChecked():
                preview_lines.append(f"  • {base_filename}_cdc_analysis.json")
            if self.format_txt.isChecked():
                preview_lines.append(f"  • {base_filename}_cdc_report.txt")
            preview_lines.append("")
        
        if self.export_visualization.isChecked():
            preview_lines.append("Visualization:")
            if self.format_html.isChecked():
                preview_lines.append(f"  • {base_filename}_viewer.html")
            if self.format_pymol.isChecked():
                preview_lines.append(f"  • {base_filename}_structure.pml")
            if self.format_png.isChecked():
                preview_lines.append(f"  • {base_filename}_plots/")
            preview_lines.append("")
        
        if self.export_summary.isChecked():
            preview_lines.append("Summary:")
            preview_lines.append(f"  • {base_filename}_summary.txt")
            preview_lines.append("")
        
        if self.compress_output.isChecked():
            preview_lines.append(f"📦 All files will be compressed into: {base_filename}.zip")
        
        self.preview_text.setText("\n".join(preview_lines))
    
    def get_export_settings(self) -> Dict[str, Any]:
        """Get the current export settings"""
        return {
            'output_dir': self.output_dir_edit.text(),
            'base_filename': self.base_filename_edit.text(),
            'export_site_energies': self.export_site_energies.isChecked() and self.export_site_energies.isEnabled(),
            'export_cdc': self.export_cdc.isChecked() and self.export_cdc.isEnabled(),
            'export_visualization': self.export_visualization.isChecked(),
            'export_summary': self.export_summary.isChecked(),
            'format_csv': self.format_csv.isChecked(),
            'format_json': self.format_json.isChecked(),
            'format_txt': self.format_txt.isChecked(),
            'format_png': self.format_png.isChecked(),
            'format_html': self.format_html.isChecked(),
            'format_pymol': self.format_pymol.isChecked(),
            'cdc_cutoff': self.cdc_cutoff_spin.value(),
            'decimal_places': self.decimal_places_spin.value(),
            'include_metadata': self.include_metadata.isChecked(),
            'compress_output': self.compress_output.isChecked()
        }
    
    def accept(self):
        """Validate and accept the dialog"""
        output_dir = self.output_dir_edit.text().strip()
        base_filename = self.base_filename_edit.text().strip()
        
        if not output_dir:
            self.preview_text.setText("❌ Error: Please specify an output directory.")
            return
        
        if not base_filename:
            self.preview_text.setText("❌ Error: Please specify a base filename.")
            return
        
        # Check if at least one content type is selected
        if not any([
            self.export_site_energies.isChecked() and self.export_site_energies.isEnabled(),
            self.export_cdc.isChecked() and self.export_cdc.isEnabled(),
            self.export_visualization.isChecked(),
            self.export_summary.isChecked()
        ]):
            self.preview_text.setText("❌ Error: Please select at least one content type to export.")
            return
        
        # Check if at least one format is selected
        if not any([
            self.format_csv.isChecked(),
            self.format_json.isChecked(),
            self.format_txt.isChecked(),
            self.format_png.isChecked(),
            self.format_html.isChecked(),
            self.format_pymol.isChecked()
        ]):
            self.preview_text.setText("❌ Error: Please select at least one export format.")
            return
        
        super().accept()
