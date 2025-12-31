"""
3D Protein Viewer Widget using py3Dmol
"""

import json
import tempfile
import os
from typing import Dict, Any, Optional, List
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton, 
                            QLabel, QComboBox, QSlider, QCheckBox, QGroupBox,
                            QSplitter, QTextEdit, QListWidget, QListWidgetItem,
                            QFrame, QScrollArea)
from PyQt5.QtCore import Qt, pyqtSignal, QUrl, QTimer
from PyQt5.QtGui import QFont
# Try to import QtWebEngineWidgets, fall back to simple viewer if not available
try:
    from PyQt5.QtWebEngineWidgets import QWebEngineView
    WEBENGINE_AVAILABLE = True
except ImportError:
    WEBENGINE_AVAILABLE = False
    QWebEngineView = None

# Import Alprotein components
from ...core.protein_structure import ProteinStructure
from ...core.pigment_system import PigmentSystem
from ...visualization.protein_3d import Protein3DRenderer


class ControlPanel(QWidget):
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


class Viewer3D(QWidget):
    """3D viewer - uses web engine if available, otherwise shows placeholder"""
    
    def __init__(self):
        super().__init__()
        self.renderer = None
        self.temp_files = []  # Track temporary files for cleanup
        
        # Set minimum size
        self.setMinimumSize(600, 400)
        
        # Setup layout
        layout = QVBoxLayout(self)
        
        if WEBENGINE_AVAILABLE:
            # Use web engine for 3D visualization
            self.web_view = QWebEngineView()
            layout.addWidget(self.web_view)
            self.load_empty_view()
        else:
            # Fallback to simple text display
            self.setup_fallback_view(layout)
    
    def setup_fallback_view(self, layout):
        """Setup fallback view when web engine is not available"""
        # Info label
        info_label = QLabel("""🧬 3D Visualization
        
The interactive 3D viewer requires QtWebEngineWidgets.
Install it with: pip install PyQtWebEngine
        
Alternatively, use the export function to generate:
• PyMOL scripts (.pml files)
• Static images
• HTML files for external viewing""")
        info_label.setStyleSheet("""
            QLabel {
                background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                          stop:0 #f7f7f7, stop:1 #efefef);
                color: #111111;
                padding: 20px;
                border-radius: 2px;
                font-size: 14px;
                text-align: center;
            }
        """)
        info_label.setAlignment(Qt.AlignCenter)
        info_label.setWordWrap(True)
        layout.addWidget(info_label)
        
        # Structure info area
        self.structure_info = QTextEdit()
        self.structure_info.setReadOnly(True)
        self.structure_info.setMaximumHeight(200)
        self.structure_info.setStyleSheet("""
            QTextEdit {
                background-color: #f3f3f3;
                border: 1px solid #e1e1e1;
                border-radius: 2px;
                padding: 10px;
                font-family: "Courier New", monospace;
                font-size: 12px;
            }
        """)
        self.structure_info.setPlainText("No structure loaded")
        layout.addWidget(self.structure_info)
        
        # Export button
        export_button = QPushButton("📁 Export for External Visualization")
        export_button.setStyleSheet("""
            QPushButton {
                background-color: #111111;
                color: #ffffff;
                border: 1px solid #111111;
                padding: 10px 20px;
                border-radius: 2px;
                font-size: 14px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #000000;
            }
        """)
        export_button.clicked.connect(self.export_visualization)
        layout.addWidget(export_button)
    
    def export_visualization(self):
        """Export visualization for external viewing"""
        if hasattr(self, 'renderer') and self.renderer:
            from PyQt5.QtWidgets import QFileDialog, QMessageBox
            
            output_dir = QFileDialog.getExistingDirectory(
                self, "Select Export Directory"
            )
            
            if output_dir:
                try:
                    exported_files = self.renderer.export_visualization_data(output_dir)
                    # QMessageBox.information( #TODO: it need to be just a message print
                    #     self, "Export Complete",
                    #     f"Visualization files exported to:\n{output_dir}\n\n" +
                    #     "\n".join([f"• {os.path.basename(f)}" for f in exported_files.values()])
                    # )
                except Exception as e:
                    QMessageBox.critical(self, "Export Error", f"Failed to export: {str(e)}")
        else:
            from PyQt5.QtWidgets import QMessageBox
            QMessageBox.warning(self, "No Structure", "No structure loaded to export.")
    
    def load_empty_view(self):
        """Load an empty 3D view"""
        if not WEBENGINE_AVAILABLE:
            return
            
        html_content = """
        <!DOCTYPE html>
        <html>
        <head>
            <title>Alprotein 3D Viewer</title>
            <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
            <style>
                body { margin: 0; padding: 0; background: #f7f7f7; font-family: Arial, sans-serif; }
                #viewer { width: 100%; height: 100vh; background: linear-gradient(45deg, #f7f7f7 0%, #efefef 100%); }
                #info { position: absolute; top: 20px; left: 20px; color: #111111; font-size: 18px; z-index: 1000; }
            </style>
        </head>
        <body>
            <div id="info">🧬 Load a PDB file to begin</div>
            <div id="viewer"></div>
            <script>
                let viewer = $3Dmol.createViewer('viewer', {
                    defaultcolors: $3Dmol.elementColors.rasmol
                });
                viewer.setBackgroundColor('#f7f7f7');
                viewer.render();
            </script>
        </body>
        </html>
        """
        self.web_view.setHtml(html_content)
    
    def load_structure(self, structure: ProteinStructure, pigment_system: PigmentSystem):
        """Load protein structure into the 3D viewer"""
        try:
            # Create renderer
            self.renderer = Protein3DRenderer(structure, pigment_system)
            
            if WEBENGINE_AVAILABLE:
                # Generate HTML with embedded structure
                html_content = self.renderer.generate_html_view()
                
                # Load the HTML
                self.web_view.setHtml(html_content)

                # Once the page finishes loading, reset camera so the model is visible
                def _on_loaded(ok):
                    self.reset_view()
                    try:
                        self.web_view.loadFinished.disconnect(_on_loaded)
                    except Exception:
                        pass

                self.web_view.loadFinished.connect(_on_loaded)
            else:
                # Update fallback view with structure info
                self.update_fallback_structure_info(structure, pigment_system)
            
        except Exception as e:
            print(f"Error loading structure in 3D viewer: {e}")
            self.load_error_view(str(e))

    def reset_view(self):
        """Reset the 3D view/camera to ensure the model is visible"""
        if not WEBENGINE_AVAILABLE or not hasattr(self, "web_view"):
            return

        if not self.web_view or not self.web_view.page():
            return

        def _run_reset():
            js_code = """
                if (typeof resetView === 'function') {
                    resetView();
                } else if (window.alproteinViewer) {
                    alproteinViewer.zoomTo();
                    alproteinViewer.render();
                }
            """
            self.web_view.page().runJavaScript(js_code)

        # Delay slightly so the HTML/JS has time to mount before we call into it
        QTimer.singleShot(100, _run_reset)
    
    def update_fallback_structure_info(self, structure: ProteinStructure, pigment_system: PigmentSystem):
        """Update structure info in fallback mode"""
        if hasattr(self, 'structure_info'):
            info_lines = [
                f"Structure: {structure.name}",
                f"Pigments: {len(pigment_system.pigments)}",
                f"Chains: {len(pigment_system.get_chains())}",
                "",
                "Pigment Details:"
            ]
            
            for pigment_id, pigment in pigment_system.pigments.items():
                info_lines.append(f"  • {pigment_id}: {pigment.get_resname()}")
            
            info_lines.extend([
                "",
                "Use the export button below to generate files for",
                "external visualization with PyMOL or other tools."
            ])
            
            self.structure_info.setPlainText("\n".join(info_lines))
    
    def load_error_view(self, error_message: str):
        """Load error view"""
        if WEBENGINE_AVAILABLE:
            html_content = f"""
            <!DOCTYPE html>
            <html>
            <head>
                <title>Alprotein 3D Viewer - Error</title>
                <style>
                    body {{ margin: 0; padding: 0; background: #f7f7f7; font-family: Arial, sans-serif; }}
                    #viewer {{ width: 100%; height: 100vh; background: #efefef; display: flex; align-items: center; justify-content: center; }}
                    #error {{ color: #111111; text-align: center; max-width: 600px; padding: 20px; }}
                </style>
            </head>
            <body>
                <div id="viewer">
                    <div id="error">
                        <h2>❌ Visualization Error</h2>
                        <p>{error_message}</p>
                        <p>Please check the structure file and try again.</p>
                    </div>
                </div>
            </body>
            </html>
            """
            self.web_view.setHtml(html_content)
        else:
            if hasattr(self, 'structure_info'):
                self.structure_info.setPlainText(f"Error loading structure:\n{error_message}")
    
    def update_view_settings(self, settings: Dict[str, Any]):
        """Update 3D view settings"""
        if not self.renderer or not WEBENGINE_AVAILABLE:
            return
        
        try:
            # Handle specific actions
            if 'action' in settings:
                if settings['action'] == 'focus_pigment':
                    pigment_id = settings.get('pigment_id', '')
                    js_code = f"""
                    viewer.zoomTo({{resi: '{pigment_id}'}});
                    viewer.render();
                    """
                    self.web_view.page().runJavaScript(js_code)
                return
            
            # Update visual styles
            js_commands = []
            
            # Update protein style
            if settings.get('show_protein', True):
                protein_style = settings.get('protein_style', 'cartoon').lower()
                js_commands.append(f"viewer.setStyle({{not: {{resn: ['CLA', 'CHL', 'BCL']}}}}, {{{protein_style}: {{}}}});")
            else:
                js_commands.append("viewer.setStyle({not: {resn: ['CLA', 'CHL', 'BCL']}}, {});")
            
            # Update pigment style
            if settings.get('show_pigments', True):
                pigment_style = settings.get('pigment_style', 'stick').lower().replace(' & ', '')
                js_commands.append(f"viewer.setStyle({{resn: ['CLA', 'CHL', 'BCL']}}, {{{pigment_style}: {{}}}});")
            else:
                js_commands.append("viewer.setStyle({resn: ['CLA', 'CHL', 'BCL']}, {});")
            
            # Apply changes
            js_commands.append("viewer.render();")
            
            # Execute JavaScript
            self.web_view.page().runJavaScript("\n".join(js_commands))
            
        except Exception as e:
            print(f"Error updating 3D view: {e}")
    
    def update_energy_visualization(self, site_energies: Dict[str, float]):
        """Update visualization with site energy data"""
        if not self.renderer:
            return

        if WEBENGINE_AVAILABLE:
            try:
                # Generate energy color mapping
                energy_colors = self.renderer.generate_energy_colors(site_energies)

                # Apply energy colors via JavaScript
                js_commands = []
                for pigment_id, color_hex in energy_colors.items():
                    js_commands.append(f"""
                    viewer.setStyle({{resi: '{pigment_id}'}}, {{stick: {{color: '{color_hex}'}}}});
                    """)

                js_commands.append("viewer.render();")
                self.web_view.page().runJavaScript("\n".join(js_commands))

            except Exception as e:
                print(f"Error updating energy visualization: {e}")
        else:
            # Update fallback view with energy information
            if hasattr(self, 'structure_info'):
                current_text = self.structure_info.toPlainText()
                energy_text = "\n\nSite Energies:\n"
                for pigment_id, energy in site_energies.items():
                    energy_text += f"  • {pigment_id}: {energy:.1f} cm⁻¹\n"

                self.structure_info.setPlainText(current_text + energy_text)

    def update_domain_visualization(self, domains_dict: Dict[int, List[int]]):
        """Update visualization with domain-based coloring"""
        if not self.renderer:
            return

        if WEBENGINE_AVAILABLE:
            try:
                # Generate domain color mapping
                domain_colors = self.renderer.generate_domain_colors(domains_dict)

                # Get list of pigment IDs to extract residue numbers
                pigment_ids = list(self.renderer.pigment_system.pigments.keys())

                # Apply domain colors via JavaScript
                js_commands = []
                for pigment_id, color_hex in domain_colors.items():
                    # Extract residue number from pigment ID (format: Chain_ResName_ResID)
                    residue_id = pigment_id.split('_')[2]
                    js_commands.append(f"""
                    viewer.setStyle({{resi: {residue_id}}}, {{stick: {{radius: 0.22, color: '{color_hex}'}}, sphere: {{radius: 0.28, opacity: 0.7, color: '{color_hex}'}}}});
                    """)

                js_commands.append("viewer.render();")
                self.web_view.page().runJavaScript("\n".join(js_commands))

            except Exception as e:
                print(f"Error updating domain visualization: {e}")
        else:
            # Update fallback view with domain information
            if hasattr(self, 'structure_info'):
                current_text = self.structure_info.toPlainText()
                domain_text = "\n\nDomains:\n"
                for domain_id, pigment_indices in domains_dict.items():
                    domain_text += f"  • Domain {domain_id + 1}: {len(pigment_indices)} pigments\n"

                self.structure_info.setPlainText(current_text + domain_text)
    
    def cleanup(self):
        """Clean up temporary files"""
        for temp_file in self.temp_files:
            try:
                os.unlink(temp_file)
            except:
                pass
        self.temp_files.clear()


class ProteinViewer(QWidget):
    """
    Main 3D protein viewer widget combining controls and visualization
    """
    pigment_highlighted = pyqtSignal(str)  # pigment_id
    
    def __init__(self):
        super().__init__()
        self.structure = None
        self.pigment_system = None
        self.site_energies = None
        self.contributions = None
        
        self.setup_ui()
        # No direct signal connections here anymore, handled by MainWindow
    
    def setup_ui(self):
        """Setup the UI layout"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # 3D viewer
        self.viewer_3d = Viewer3D()
        layout.addWidget(self.viewer_3d)
    
    # Removed connect_signals as it's no longer needed here
    
    def load_structure(self, structure: ProteinStructure, pigment_system: PigmentSystem):
        """Load protein structure and pigment system"""
        self.structure = structure
        self.pigment_system = pigment_system
        
        # Load into 3D viewer
        self.viewer_3d.load_structure(structure, pigment_system)

        # Auto-reset view so the freshly loaded structure is visible without user action
        self.reset_view()
    
    def update_energy_visualization(self, site_energies: Dict[str, float]):
        """Update visualization with site energy results"""
        self.site_energies = site_energies
        self.viewer_3d.update_energy_visualization(site_energies)
    
    def update_contribution_visualization(self, contributions: Dict[str, Dict[str, Any]]):
        """Update visualization with CDC contribution results"""
        self.contributions = contributions
        # TODO: Implement contribution-based coloring

    def update_domain_visualization(self, domains_dict: Dict[int, List[int]]):
        """Update visualization with domain-based coloring"""
        self.viewer_3d.update_domain_visualization(domains_dict)

    def highlight_pigment(self, pigment_id: str):
        """Highlight a specific pigment"""
        settings = {'action': 'focus_pigment', 'pigment_id': pigment_id}
        self.viewer_3d.update_view_settings(settings)

    def reset_view(self):
        """Reset the underlying 3D view if available"""
        if self.viewer_3d:
            self.viewer_3d.reset_view()
    
    def cleanup(self):
        """Clean up resources"""
        if self.viewer_3d:
            self.viewer_3d.cleanup()
