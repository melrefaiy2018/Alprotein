"""
Data Table Widget for Scientific Workbench

Professional data table with filtering, sorting, and export capabilities.
"""

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTableWidget, QTableWidgetItem, QHeaderView, QComboBox,
    QLineEdit, QFrame, QSizePolicy, QMenu, QAction
)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont, QColor, QBrush, QCursor
from typing import Dict, List, Any
import numpy as np


class DataTableWidget(QWidget):
    """
    Professional data table for displaying calculation results
    """

    row_selected = pyqtSignal(str)  # pigment_id
    export_requested = pyqtSignal()
    plot_requested = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.data = {}
        self.setup_ui()

    def setup_ui(self):
        """Setup the UI"""
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(16, 16, 16, 16)
        main_layout.setSpacing(16)

        # Header with title and controls
        header_layout = QHBoxLayout()

        title = QLabel("DATA TABLE")
        title.setStyleSheet("""
            QLabel {
                font-size: 13px;
                font-weight: bold;
                color: #111111;
            }
        """)
        header_layout.addWidget(title)

        header_layout.addStretch()

        # Filter dropdown
        self.filter_combo = QComboBox()
        self.filter_combo.addItems(["All Data", "Site Energies Only", "High Energy", "Low Energy"])
        self.filter_combo.setStyleSheet("""
            QComboBox {
                padding: 6px;
                border: 1px solid #e1e1e1;
                border-radius: 2px;
                background-color: #ffffff;
                color: #111111;
                min-width: 150px;
            }
        """)
        self.filter_combo.currentTextChanged.connect(self.apply_filter)
        header_layout.addWidget(self.filter_combo)

        # Search box
        self.search_box = QLineEdit()
        self.search_box.setPlaceholderText("Search...")
        self.search_box.setMaximumWidth(200)
        self.search_box.setStyleSheet("""
            QLineEdit {
                padding: 6px;
                border: 1px solid #e1e1e1;
                border-radius: 2px;
                background-color: #ffffff;
                color: #111111;
            }
            QLineEdit:focus {
                border: 2px solid #111111;
            }
        """)
        self.search_box.textChanged.connect(self.apply_filter)
        header_layout.addWidget(self.search_box)

        main_layout.addLayout(header_layout)

        # Table
        self.table = QTableWidget()
        self.table.setColumnCount(5)
        self.table.setHorizontalHeaderLabels([
            "ID", "Site Energy (cm⁻¹)", "Shift (cm⁻¹)",
            "Coupling Range", "Contributors"
        ])

        # Table styling
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setSelectionMode(QTableWidget.SingleSelection)
        self.table.setSortingEnabled(True)
        self.table.setStyleSheet("""
            QTableWidget {
                background-color: #ffffff;
                alternate-background-color: #f7f7f7;
                gridline-color: #e6e6e6;
                border: 1px solid #e1e1e1;
                border-radius: 2px;
                color: #111111;
            }
            QHeaderView::section {
                background-color: #111111;
                color: #ffffff;
                padding: 8px;
                border: none;
                font-weight: bold;
                font-size: 11px;
            }
            QTableWidget::item {
                padding: 6px;
                color: #111111;
            }
            QTableWidget::item:selected {
                background-color: #111111;
                color: #ffffff;
            }
        """)

        # Set column widths
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(3, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(4, QHeaderView.Stretch)

        # Connect selection signal
        self.table.itemSelectionChanged.connect(self.on_selection_changed)

        # Context menu
        self.table.setContextMenuPolicy(Qt.CustomContextMenu)
        self.table.customContextMenuRequested.connect(self.show_context_menu)

        main_layout.addWidget(self.table)

        # Action buttons
        button_layout = QHBoxLayout()

        self.copy_btn = QPushButton("📋 Copy Selected")
        self.copy_btn.setStyleSheet(self._button_style("#111111"))
        self.copy_btn.clicked.connect(self.copy_selected)
        button_layout.addWidget(self.copy_btn)

        self.export_csv_btn = QPushButton("💾 Export CSV")
        self.export_csv_btn.setStyleSheet(self._button_style("#111111"))
        self.export_csv_btn.clicked.connect(self.export_requested.emit)
        button_layout.addWidget(self.export_csv_btn)

        self.export_json_btn = QPushButton("📄 Export JSON")
        self.export_json_btn.setStyleSheet(self._button_style("#111111"))
        self.export_json_btn.clicked.connect(self.export_json)
        button_layout.addWidget(self.export_json_btn)

        self.plot_btn = QPushButton("📊 Plot")
        self.plot_btn.setStyleSheet(self._button_style("#111111"))
        self.plot_btn.clicked.connect(self.plot_requested.emit)
        button_layout.addWidget(self.plot_btn)

        button_layout.addStretch()

        main_layout.addLayout(button_layout)

        # Status label
        self.status_label = QLabel("No data loaded")
        self.status_label.setStyleSheet("color: #6b6b6b; font-style: italic; font-size: 11px;")
        main_layout.addWidget(self.status_label)

    def _button_style(self, color):
        """Generate button style with given color"""
        return f"""
            QPushButton {{
                background-color: {color};
                color: white;
                border: 1px solid {color};
                padding: 8px 16px;
                border-radius: 2px;
                font-weight: bold;
            }}
            QPushButton:hover {{
                background-color: {self._darken_color(color)};
            }}
            QPushButton:pressed {{
                background-color: {self._darken_color(color, 0.3)};
            }}
        """

    def _darken_color(self, hex_color, factor=0.2):
        """Darken a hex color"""
        # Simple darkening by reducing RGB values
        r = int(hex_color[1:3], 16)
        g = int(hex_color[3:5], 16)
        b = int(hex_color[5:7], 16)

        r = int(r * (1 - factor))
        g = int(g * (1 - factor))
        b = int(b * (1 - factor))

        return f"#{r:02x}{g:02x}{b:02x}"

    def update_data(self, site_energies: Dict[str, float], vacuum_energies: Dict[str, float],
                    couplings: Dict[str, Any] = None, contributions: Dict[str, Any] = None):
        """Update table with calculation results"""
        self.data = {
            'site_energies': site_energies,
            'vacuum_energies': vacuum_energies,
            'couplings': couplings or {},
            'contributions': contributions or {}
        }

        self.populate_table()
        self.status_label.setText(f"Showing {len(site_energies)} pigments")
        self.status_label.setStyleSheet("color: #1f7a44; font-weight: bold; font-size: 11px;")

    def populate_table(self):
        """Populate table with data"""
        self.table.setRowCount(0)
        self.table.setSortingEnabled(False)

        site_energies = self.data.get('site_energies', {})
        vacuum_energies = self.data.get('vacuum_energies', {})
        couplings = self.data.get('couplings', {})
        contributions = self.data.get('contributions', {})

        for i, (pigment_id, site_energy) in enumerate(site_energies.items()):
            self.table.insertRow(i)

            # ID
            id_item = QTableWidgetItem(pigment_id)
            id_item.setFont(QFont("Monospace", 10, QFont.Bold))
            self.table.setItem(i, 0, id_item)

            # Site Energy
            se_item = QTableWidgetItem(f"{site_energy:.1f}")
            se_item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            se_item.setData(Qt.UserRole, site_energy)
            self.table.setItem(i, 1, se_item)

            # Shift
            vacuum_e = vacuum_energies.get(pigment_id, 0)
            shift = site_energy - vacuum_e
            shift_item = QTableWidgetItem(f"{shift:+.1f}")
            shift_item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            shift_item.setData(Qt.UserRole, shift)

            # Color code by shift
            if shift > 0:
                shift_item.setBackground(QBrush(QColor(255, 230, 230)))  # Light red
            elif shift < 0:
                shift_item.setBackground(QBrush(QColor(230, 255, 230)))  # Light green

            self.table.setItem(i, 2, shift_item)

            # Coupling range (placeholder for now)
            coupling_item = QTableWidgetItem("-")
            coupling_item.setTextAlignment(Qt.AlignCenter)
            self.table.setItem(i, 3, coupling_item)

            # Contributors (placeholder)
            contrib_text = "Protein"
            if pigment_id in contributions:
                contrib_data = contributions[pigment_id]
                if isinstance(contrib_data, dict):
                    # Get top contributors
                    top_contrib = sorted(contrib_data.items(),
                                       key=lambda x: abs(x[1]) if isinstance(x[1], (int, float)) else 0,
                                       reverse=True)[:2]
                    if top_contrib:
                        contrib_text = ", ".join([k for k, v in top_contrib])

            contrib_item = QTableWidgetItem(contrib_text)
            self.table.setItem(i, 4, contrib_item)

        self.table.setSortingEnabled(True)

    def apply_filter(self):
        """Apply filter to table"""
        filter_text = self.filter_combo.currentText()
        search_text = self.search_box.text().lower()

        for row in range(self.table.rowCount()):
            show_row = True

            # Apply search filter
            if search_text:
                row_text = ""
                for col in range(self.table.columnCount()):
                    item = self.table.item(row, col)
                    if item:
                        row_text += item.text().lower() + " "

                if search_text not in row_text:
                    show_row = False

            # Apply category filter
            if show_row and filter_text != "All Data":
                shift_item = self.table.item(row, 2)
                if shift_item:
                    shift = shift_item.data(Qt.UserRole)
                    if filter_text == "High Energy" and shift <= 0:
                        show_row = False
                    elif filter_text == "Low Energy" and shift >= 0:
                        show_row = False

            self.table.setRowHidden(row, not show_row)

    def on_selection_changed(self):
        """Handle row selection"""
        selected = self.table.selectedItems()
        if selected:
            row = selected[0].row()
            pigment_id = self.table.item(row, 0).text()
            self.row_selected.emit(pigment_id)

    def show_context_menu(self, position):
        """Show context menu"""
        menu = QMenu()

        copy_action = QAction("Copy Row", self)
        copy_action.triggered.connect(self.copy_selected)
        menu.addAction(copy_action)

        export_action = QAction("Export Selected", self)
        export_action.triggered.connect(self.export_selected)
        menu.addAction(export_action)

        menu.exec_(self.table.viewport().mapToGlobal(position))

    def copy_selected(self):
        """Copy selected row to clipboard"""
        from PyQt5.QtWidgets import QApplication

        selected = self.table.selectedItems()
        if selected:
            row_data = []
            row = selected[0].row()
            for col in range(self.table.columnCount()):
                item = self.table.item(row, col)
                if item:
                    row_data.append(item.text())

            QApplication.clipboard().setText("\t".join(row_data))
            self.status_label.setText("Row copied to clipboard")

    def export_selected(self):
        """Export selected row"""
        # Could implement specific export for selected row
        self.export_requested.emit()

    def export_json(self):
        """Export data as JSON"""
        from PyQt5.QtWidgets import QFileDialog, QMessageBox
        import json
        import sys

        if not self.data:
            return

        options = QFileDialog.DontUseNativeDialog if sys.platform == 'darwin' else QFileDialog.Options()

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export as JSON",
            "results.json",
            "JSON Files (*.json);;All Files (*)",
            options=options
        )

        if file_path:
            try:
                # Convert numpy types to Python types
                export_data = {}
                for key, value in self.data.items():
                    if isinstance(value, dict):
                        export_data[key] = {
                            k: float(v) if isinstance(v, (np.floating, np.integer)) else v
                            for k, v in value.items()
                        }
                    else:
                        export_data[key] = value

                with open(file_path, 'w') as f:
                    json.dump(export_data, f, indent=2)

                QMessageBox.information(self, "Export Successful",
                                      f"Data exported to:\n{file_path}")
            except Exception as e:
                QMessageBox.critical(self, "Export Error",
                                   f"Failed to export:\n{str(e)}")

    def clear_data(self):
        """Clear all data"""
        self.table.setRowCount(0)
        self.data = {}
        self.status_label.setText("No data loaded")
        self.status_label.setStyleSheet("color: #6b6b6b; font-style: italic; font-size: 11px;")
