"""
File Loader Widget for PDB files
"""

import os
import sys
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton, 
                            QLabel, QGroupBox, QFileDialog, QTextEdit, QFrame)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont, QPixmap, QPainter, QPen, QDragEnterEvent, QDropEvent


class DropZone(QFrame):
    """
    Drag and drop zone for PDB files
    """
    file_dropped = pyqtSignal(str)
    
    def __init__(self):
        super().__init__()
        self.setAcceptDrops(True)
        self.setMinimumHeight(80)
        self.setFrameStyle(QFrame.StyledPanel)
        self.setStyleSheet("""
            DropZone {
                border: 1px dashed #e1e1e1;
                border-radius: 2px;
                background-color: #f7f7f7;
            }
            DropZone:hover {
                border-color: #111111;
                background-color: #efefef;
            }
        """)
        
        layout = QVBoxLayout(self)
        layout.setAlignment(Qt.AlignCenter)
        
        # Icon and text
        self.label = QLabel("📁\n\nDrag & Drop PDB File Here\nor click Browse below")
        self.label.setAlignment(Qt.AlignCenter)
        self.label.setStyleSheet("""
            QLabel {
                color: #6b6b6b;
                font-size: 14px;
                border: none;
                background: transparent;
            }
        """)
        layout.addWidget(self.label)
    
    def dragEnterEvent(self, event: QDragEnterEvent):
        """Handle drag enter event"""
        if event.mimeData().hasUrls():
            urls = event.mimeData().urls()
            if len(urls) == 1 and urls[0].toLocalFile().endswith('.pdb'):
                event.accept()
                self.setStyleSheet("""
                    DropZone {
                        border: 2px dashed #1f7a44;
                        border-radius: 2px;
                        background-color: #eaf3ed;
                    }
                """)
                return
        event.ignore()
    
    def dragLeaveEvent(self, event):
        """Handle drag leave event"""
        self.setStyleSheet("""
            DropZone {
                border: 1px dashed #e1e1e1;
                border-radius: 2px;
                background-color: #f7f7f7;
            }
            DropZone:hover {
                border-color: #111111;
                background-color: #efefef;
            }
        """)
    
    def dropEvent(self, event: QDropEvent):
        """Handle drop event"""
        urls = event.mimeData().urls()
        if urls:
            file_path = urls[0].toLocalFile()
            if file_path.endswith('.pdb'):
                self.file_dropped.emit(file_path)
                event.accept()
            else:
                event.ignore()
        
        # Reset styling
        self.dragLeaveEvent(event)


class FileLoaderWidget(QWidget):
    """
    Widget for loading PDB files with drag-and-drop and file browser
    """
    file_selected = pyqtSignal(str)  # file_path
    
    def __init__(self):
        super().__init__()
        self.current_file = None
        self.setup_ui()
    
    def setup_ui(self):
        """Setup the UI components"""
        layout = QVBoxLayout(self)
        
        # Set widget styling to ensure black text
        self.setStyleSheet("""
            QWidget {
                color: black;
            }
            QLabel {
                color: black;
            }
            QGroupBox {
                color: black;
                font-weight: bold;
            }
        """)
        
        # Create group box
        group_box = QGroupBox("📁 Load PDB Structure")
        group_layout = QVBoxLayout(group_box)
        
        # Drop zone
        self.drop_zone = DropZone()
        self.drop_zone.file_dropped.connect(self.load_file)
        group_layout.addWidget(self.drop_zone)
        
        # Browse button
        browse_layout = QHBoxLayout()
        self.browse_button = QPushButton("Browse Files...")
        self.browse_button.clicked.connect(self.browse_files)
        self.browse_button.setStyleSheet("""
            QPushButton {
                background-color: #111111;
                color: #ffffff;
                border: 1px solid #111111;
                padding: 8px 16px;
                border-radius: 2px;
                font-size: 14px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #000000;
            }
            QPushButton:pressed {
                background-color: #000000;
            }
        """)
        browse_layout.addWidget(self.browse_button)
        browse_layout.addStretch()
        group_layout.addLayout(browse_layout)
        
        # File info area
        self.file_info = QTextEdit()
        self.file_info.setMaximumHeight(60)
        self.file_info.setReadOnly(True)
        self.file_info.setVisible(False)
        self.file_info.setStyleSheet("""
            QTextEdit {
                color: #111111;
                background-color: #f3f3f3;
                border: 1px solid #e1e1e1;
                border-radius: 2px;
                padding: 8px;
                font-family: "Courier New", monospace;
                font-size: 12px;
            }
        """)
        group_layout.addWidget(self.file_info)
        
        layout.addWidget(group_box)
        layout.addStretch()
    
    def browse_files(self):
        """Open file browser dialog"""
        # Use non-native dialog on macOS to avoid warnings
        options = QFileDialog.DontUseNativeDialog if sys.platform == 'darwin' else QFileDialog.Options()
        
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select PDB File",
            "",
            "PDB Files (*.pdb);;All Files (*)",
            options=options
        )
        
        if file_path:
            self.load_file(file_path)
    
    def load_file(self, file_path: str):
        """Load the selected file"""
        if not os.path.exists(file_path):
            return
        
        self.current_file = file_path
        self.file_selected.emit(file_path)
        
        # Update UI to show loading state
        self.drop_zone.label.setText("📁\n\nLoading...")
        self.browse_button.setEnabled(False)
    
    def set_loaded_file(self, file_path: str, pigment_count: int):
        """Update UI after successful file loading"""
        self.current_file = file_path
        filename = os.path.basename(file_path)
        file_size = os.path.getsize(file_path)
        
        # Update drop zone
        self.drop_zone.label.setText(f"✅\n\n{filename}\nLoaded Successfully")
        self.drop_zone.setStyleSheet("""
            DropZone {
                border: 2px solid #1f7a44;
                border-radius: 2px;
                background-color: #eaf3ed;
            }
        """)
        
        # Show file info
        info_text = f"""File: {filename}
Size: {file_size / 1024:.1f} KB
Pigments: {pigment_count}
Path: {file_path}"""
        
        self.file_info.setText(info_text)
        self.file_info.setVisible(True)
        
        # Re-enable browse button
        self.browse_button.setEnabled(True)
        self.browse_button.setText("Load Different File...")
    
    def set_loading_error(self, error_message: str):
        """Update UI after loading error"""
        self.drop_zone.label.setText("❌\n\nLoading Failed\nClick Browse to try again")
        self.drop_zone.setStyleSheet("""
            DropZone {
                border: 2px dashed #8b2b2b;
                border-radius: 2px;
                background-color: #f6eeee;
            }
        """)
        
        # Show error info
        self.file_info.setText(f"Error: {error_message}")
        self.file_info.setVisible(True)
        
        # Re-enable browse button
        self.browse_button.setEnabled(True)
        self.browse_button.setText("Browse Files...")
    
    def clear(self):
        """Clear the current file and reset UI"""
        self.current_file = None
        
        # Reset drop zone
        self.drop_zone.label.setText("📁\n\nDrag & Drop PDB File Here\nor click Browse below")
        self.drop_zone.setStyleSheet("""
            DropZone {
                border: 1px dashed #e1e1e1;
                border-radius: 2px;
                background-color: #f7f7f7;
            }
            DropZone:hover {
                border-color: #111111;
                background-color: #efefef;
            }
        """)
        
        # Hide file info
        self.file_info.setVisible(False)
        
        # Reset browse button
        self.browse_button.setText("Browse Files...")
        self.browse_button.setEnabled(True)
