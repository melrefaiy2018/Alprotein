from typing import Dict, Any, List, Optional
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTabWidget,
    QTableWidget, QTableWidgetItem, QHeaderView,
    QPushButton, QLabel, QTextEdit, QGroupBox,
    QSplitter, QScrollArea, QFrame, QSizePolicy,
    QComboBox, QSpinBox, QCheckBox
)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont, QColor, QBrush

class ContributionAnalysisWidget(QWidget):
    """Widget for displaying CDC contribution analysis"""
    
    pigment_selected = pyqtSignal(str)  # pigment_id
    
    def __init__(self):
        super().__init__()
        self.contributions = {}
        self.current_pigment = None
        self.setup_ui()
    
    def setup_ui(self):
        """Setup the UI layout"""
        layout = QVBoxLayout(self)
        
        # Controls
        controls_layout = QHBoxLayout()
        
        controls_layout.addWidget(QLabel("Pigment:"))
        self.pigment_combo = QComboBox()
        self.pigment_combo.currentTextChanged.connect(self.on_pigment_changed)
        controls_layout.addWidget(self.pigment_combo)
        
        controls_layout.addWidget(QLabel("Cutoff:"))
        self.cutoff_spin = QSpinBox()
        self.cutoff_spin.setRange(1, 1000)
        self.cutoff_spin.setValue(20)
        self.cutoff_spin.setSuffix(" cm⁻¹")
        self.cutoff_spin.valueChanged.connect(self.update_display)
        controls_layout.addWidget(self.cutoff_spin)
        
        controls_layout.addStretch()
        layout.addLayout(controls_layout)
        
        # Splitter for table and plot
        splitter = QSplitter(Qt.Vertical)
        
        # Contributions table
        self.contributions_table = QTableWidget()
        self.setup_contributions_table()
        splitter.addWidget(self.contributions_table)
        
        # Plot widget
        self.plot_widget = ContributionPlotWidget()
        splitter.addWidget(self.plot_widget)
        
        splitter.setSizes([500, 400])
        layout.addWidget(splitter)
    
    def setup_contributions_table(self):
        """Setup the contributions table"""
        self.contributions_table.setColumnCount(3)
        self.contributions_table.setHorizontalHeaderLabels([
            "Contributor", "Contribution (cm⁻¹)", "Type"
        ])
        
        header = self.contributions_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.Stretch)
        header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        
        self.contributions_table.setSortingEnabled(True)
    
    def update_data(self, contributions: Dict[str, Dict[str, Any]]):
        """Update with CDC analysis results"""
        self.contributions = contributions
        
        # Update pigment combo
        self.pigment_combo.clear()
        self.pigment_combo.addItems(sorted(contributions.keys()))
        
        if contributions:
            # Select first pigment
            first_pigment = sorted(contributions.keys())[0]
            self.pigment_combo.setCurrentText(first_pigment)
            self.on_pigment_changed(first_pigment)
    
    def on_pigment_changed(self, pigment_id: str):
        """Handle pigment selection change"""
        if not pigment_id or pigment_id not in self.contributions:
            return
        
        self.current_pigment = pigment_id
        self.pigment_selected.emit(pigment_id)
        self.update_display()
    
    def update_display(self):
        """Update the display for current pigment"""
        if not self.current_pigment or self.current_pigment not in self.contributions:
            return
        
        pigment_contribs = self.contributions[self.current_pigment]
        cutoff = self.cutoff_spin.value()
        
        # Filter contributions by cutoff
        filtered_contribs = {}
        for contributor, value in pigment_contribs.items():
            if contributor == 'vacuum':
                continue  # Skip vacuum energy in contributions display
            
            if isinstance(value, dict) and 'total_contribution' in value:
                # New format with detailed atom breakdown
                contrib_value = value['total_contribution']
            else:
                # Simple format
                contrib_value = value
            
            if abs(contrib_value) >= cutoff:
                filtered_contribs[contributor] = contrib_value
        
        # Update table
        self.update_contributions_table(filtered_contribs)
        
        # Update plot
        self.plot_widget.update_plot(filtered_contribs, self.current_pigment)
    
    def update_contributions_table(self, contributions: Dict[str, float]):
        """Update the contributions table"""
        self.contributions_table.setRowCount(0)
        
        # Sort by contribution magnitude
        sorted_contribs = sorted(contributions.items(), 
                               key=lambda x: abs(x[1]), reverse=True)
        
        for i, (contributor, value) in enumerate(sorted_contribs):
            self.contributions_table.insertRow(i)
            
            # Contributor name
            name_item = QTableWidgetItem(contributor)
            self.contributions_table.setItem(i, 0, name_item)
            
            # Contribution value
            value_item = QTableWidgetItem(f"{value:+.1f}")
            value_item.setData(Qt.UserRole, value)
            
            # Color code by value
            if value > 0:
                value_item.setBackground(QBrush(QColor(230, 230, 230)))  # Light gray
            else:
                value_item.setBackground(QBrush(QColor(210, 210, 210)))  # Darker gray
            
            self.contributions_table.setItem(i, 1, value_item)
            
            # Type
            if any(pigment_name in contributor for pigment_name in ['CLA', 'CHL', 'BCL']):
                contrib_type = "Pigment"
            elif any(ion in contributor for ion in ['NA', 'CL']):
                contrib_type = "Ion"
            else:
                contrib_type = "Residue"
            
            type_item = QTableWidgetItem(contrib_type)
            self.contributions_table.setItem(i, 2, type_item)

class ContributionPlotWidget(QWidget):
    """Widget for plotting contribution data"""
    
    def __init__(self):
        super().__init__()
        self.setup_ui()
    
    def setup_ui(self):
        """Setup the plot widget"""
        layout = QVBoxLayout(self)
        
        # Create matplotlib figure
        self.figure = Figure(figsize=(10, 4), facecolor='white')
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(self.canvas)
        
        # Plot controls
        controls_layout = QHBoxLayout()
        
        self.show_positive = QCheckBox("Show Positive")
        self.show_positive.setChecked(True)
        self.show_positive.toggled.connect(self.update_current_plot)
        controls_layout.addWidget(self.show_positive)
        
        self.show_negative = QCheckBox("Show Negative")
        self.show_negative.setChecked(True)
        self.show_negative.toggled.connect(self.update_current_plot)
        controls_layout.addWidget(self.show_negative)
        
        self.max_bars = QSpinBox()
        self.max_bars.setRange(5, 50)
        self.max_bars.setValue(20)
        self.max_bars.setPrefix("Top ")
        self.max_bars.valueChanged.connect(self.update_current_plot)
        controls_layout.addWidget(self.max_bars)
        
        controls_layout.addStretch()
        layout.addLayout(controls_layout)
        
        self.current_data = {}
        self.current_pigment = ""
    
    def update_plot(self, contributions: Dict[str, float], pigment_id: str):
        """Update the plot with new data"""
        self.current_data = contributions
        self.current_pigment = pigment_id
        self.update_current_plot()
    
    def update_current_plot(self):
        """Update the plot with current settings"""
        if not self.current_data:
            return
        
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.set_facecolor('white')
        
        # Filter data based on checkboxes
        filtered_data = {}
        for contrib, value in self.current_data.items():
            if self.show_positive.isChecked() and value > 0:
                filtered_data[contrib] = value
            elif self.show_negative.isChecked() and value < 0:
                filtered_data[contrib] = value
        
        if not filtered_data:
            ax.text(0.5, 0.5, "No data to display", 
                   ha='center', va='center', transform=ax.transAxes, color='black')
            self.canvas.draw()
            return
        
        # Sort and limit data
        sorted_data = sorted(filtered_data.items(), 
                           key=lambda x: abs(x[1]), reverse=True)
        max_items = self.max_bars.value()
        plot_data = sorted_data[:max_items]
        
        if not plot_data:
            self.canvas.draw()
            return
        
        # Prepare plot data
        names = [item[0] for item in plot_data]
        values = [item[1] for item in plot_data]
        
        # Truncate long names
        display_names = []
        for name in names:
            if len(name) > 15:
                display_names.append(name[:12] + "...")
            else:
                display_names.append(name)
        
        # Create colors
        colors = ['#6b6b6b' if v < 0 else '#111111' for v in values]
        
        # Create bar plot
        bars = ax.bar(range(len(values)), values, color=colors, alpha=0.7)
        
        # Customize plot
        ax.set_xlabel("Contributing Residues/Pigments")
        ax.set_ylabel("Energy Contribution (cm⁻¹)")
        ax.set_title(f"Site Energy Contributions for {self.current_pigment}")
        ax.set_xticks(range(len(display_names)))
        ax.set_xticklabels(display_names, rotation=45, ha='right')
        ax.grid(axis='y', alpha=0.3)
        
        # Add value labels on bars
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., 
                   height + (5 if height > 0 else -15),
                   f'{value:.0f}', 
                   ha='center', va='bottom' if height > 0 else 'top',
                   fontsize=8, color='black')
        
        # Adjust layout
        self.figure.tight_layout()
        self.canvas.draw()

class ResultsPanel(QWidget):
    """
    Redesigned results panel with a clean, modern layout using cards.
    """
    pigment_selected = pyqtSignal(str)
    export_requested = pyqtSignal()

    def __init__(self):
        super().__init__()
        self.setup_ui()
        self.connect_signals()

    def setup_ui(self):
        """Setup the main UI layout with a tabbed interface."""
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(16, 16, 16, 16)
        main_layout.setSpacing(16)

        self.tab_widget = QTabWidget()
        main_layout.addWidget(self.tab_widget)

        # Site Energies Tab
        self.site_energy_tab = QWidget()
        self.setup_site_energy_tab()
        self.tab_widget.addTab(self.site_energy_tab, "Site Energies")

        # CDC Analysis Tab
        self.cdc_tab = QWidget()
        self.setup_cdc_tab()
        self.tab_widget.addTab(self.cdc_tab, "CDC Analysis")

        # Export button
        export_layout = QHBoxLayout()
        export_layout.addStretch()
        
        self.export_button = QPushButton("📁 Export Results")
        self.export_button.setEnabled(False)
        self.export_button.clicked.connect(self.export_requested.emit)
        self.export_button.setStyleSheet("""
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
            QPushButton:disabled {
                background-color: #f3f3f3;
                color: #b0b0b0;
                border: 1px solid #e1e1e1;
            }
        """)
        export_layout.addWidget(self.export_button)
        
        main_layout.addLayout(export_layout)

    def setup_site_energy_tab(self):
        """Setup the site energies tab"""
        layout = QVBoxLayout(self.site_energy_tab)
        
        # Site energy table
        self.site_energy_table = self._create_site_energy_table({})
        layout.addWidget(self.site_energy_table)

    def setup_cdc_tab(self):
        """Setup the CDC analysis tab"""
        layout = QVBoxLayout(self.cdc_tab)
        self.cdc_widget = ContributionAnalysisWidget()
        layout.addWidget(self.cdc_widget)

    def connect_signals(self):
        """Connect signals from any interactive components."""
        pass

    def clear_results(self):
        """Clear all result cards from the layout."""
        self.site_energy_table.setRowCount(0)
        if hasattr(self, 'cdc_widget'):
            self.cdc_widget.pigment_combo.clear()
            self.cdc_widget.contributions_table.setRowCount(0)
            self.cdc_widget.plot_widget.figure.clear()
            self.cdc_widget.plot_widget.canvas.draw()

    def display_site_energies(self, results: Dict[str, Any]):
        """Display site energy results in a styled card."""
        self.clear_results()
        self._update_site_energy_table(results)

    def display_cdc_analysis(self, results: Dict[str, Any]):
        """Display CDC analysis results in a scrollable, styled card."""
        if hasattr(self, 'cdc_widget'):
            self.cdc_widget.update_data(results.get('detailed_contributions', {}))

    def _create_site_energy_table(self, results: Dict[str, Any]) -> QTableWidget:
        """Creates and populates a table for site energies."""
        table = QTableWidget()
        table.setColumnCount(4)
        table.setHorizontalHeaderLabels(["Pigment ID", "Site Energy (cm⁻¹)", "Vacuum Energy (cm⁻¹)", "Energy Shift (cm⁻¹)"])
        table.setSortingEnabled(True)
        table.setSelectionBehavior(QTableWidget.SelectRows)
        self._update_site_energy_table(results, table)
        return table

    def _update_site_energy_table(self, results: Dict[str, Any], table: QTableWidget = None):
        """Updates the site energy table with new data."""
        if table is None:
            table = self.site_energy_table

        table.setRowCount(0)
        site_energies = results.get('site_energies', {})
        vacuum_energies = results.get('vacuum_energies', {})

        for i, (pigment_id, site_energy) in enumerate(site_energies.items()):
            vacuum_energy = vacuum_energies.get(pigment_id, 0.0)
            energy_shift = site_energy - vacuum_energy

            table.insertRow(i)
            table.setItem(i, 0, QTableWidgetItem(pigment_id))
            table.setItem(i, 1, QTableWidgetItem(f"{site_energy:.1f}"))
            table.setItem(i, 2, QTableWidgetItem(f"{vacuum_energy:.1f}"))
            table.setItem(i, 3, QTableWidgetItem(f"{energy_shift:+.1f}"))

        table.resizeColumnsToContents()
