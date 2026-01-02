"""
Hamiltonian Visualization Widget for Alprotein GUI

This widget provides visualization for Hamiltonian matrix and eigenvalues:
- Interactive heatmap of coupling matrix
- Eigenvalue table with contributions
- Export capabilities
"""

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTableWidget, QTableWidgetItem, QHeaderView, QSplitter,
    QGroupBox, QCheckBox, QSpinBox, QDoubleSpinBox,
    QRadioButton, QButtonGroup, QSizePolicy, QTextEdit
)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont, QColor, QBrush

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import numpy as np
from typing import Dict, List, Tuple, Optional, Any


class HamiltonianWidget(QWidget):
    """
    Widget for visualizing Hamiltonian matrix and eigenvalues
    """

    run_hamiltonian_requested = pyqtSignal()
    export_requested = pyqtSignal()
    view_full_matrix = pyqtSignal()
    site_energy_updated = pyqtSignal(str, float)
    parameters_changed = pyqtSignal(dict)
    domains_updated = pyqtSignal(dict)  # Emitted when domain clustering is updated

    def __init__(self, parent=None):
        super().__init__(parent)
        self.hamiltonian = None
        self.hamiltonian_calculator = None  # Reference to HamiltonianCalculator
        self.eigenvalues = None
        self.eigenvectors = None
        self.pigment_labels = None
        self.site_energies = None
        self.E0a = 14900.0  # Default vacuum energy for CLA
        self._updating_site_table = False
        self.setup_ui()

    def setup_ui(self):
        """Setup the UI"""
        main_layout = QHBoxLayout(self)
        main_layout.setContentsMargins(16, 16, 16, 16)
        main_layout.setSpacing(16)

        # Left side: Site energy shift bar plot
        left_widget = self.create_site_energy_plot_widget()
        main_layout.addWidget(left_widget, stretch=3)

        # Right side: Site energies and domain information
        right_widget = self.create_info_widget()
        main_layout.addWidget(right_widget, stretch=1)

    def create_site_energy_plot_widget(self):
        """Create the site energy shift bar plot widget"""
        widget = QGroupBox("SITE ENERGY SHIFTS")
        widget.setProperty("class", "card")

        layout = QVBoxLayout(widget)
        layout.setSpacing(10)

        # Site energy shift bar plot
        self.site_energy_figure = Figure(figsize=(8, 6), dpi=100)
        self.site_energy_canvas = FigureCanvas(self.site_energy_figure)
        self.site_energy_canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        layout.addWidget(self.site_energy_canvas)

        # Plot info
        self.plot_info_label = QLabel("Site energies not calculated")
        self.plot_info_label.setProperty("class", "label")
        self.plot_info_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(self.plot_info_label)

        # Buttons
        btn_layout = QHBoxLayout()

        self.view_matrix_btn = QPushButton("📋 View Full Matrix")
        self.view_matrix_btn.clicked.connect(self.view_full_matrix.emit)
        self.view_matrix_btn.setEnabled(False)
        self.view_matrix_btn.setProperty("class", "primary")
        btn_layout.addWidget(self.view_matrix_btn)

        self.export_matrix_btn = QPushButton("💾 Export")
        self.export_matrix_btn.clicked.connect(self.export_requested.emit)
        self.export_matrix_btn.setEnabled(False)
        self.export_matrix_btn.setProperty("class", "secondary")
        btn_layout.addWidget(self.export_matrix_btn)

        btn_layout.addStretch()
        layout.addLayout(btn_layout)

        return widget

    def create_info_widget(self):
        """Create the info widget with site energies and domains"""
        widget = QGroupBox("SITE ENERGIES & DOMAINS")
        widget.setProperty("class", "card")

        layout = QVBoxLayout(widget)
        layout.setSpacing(10)

        # Hamiltonian settings
        settings_group = QGroupBox("HAMILTONIAN SETTINGS")
        settings_group.setProperty("class", "card")
        settings_layout = QVBoxLayout(settings_group)
        settings_layout.setSpacing(6)

        method_label = QLabel("Coupling Method")
        method_label.setProperty("class", "label")
        settings_layout.addWidget(method_label)

        self.method_group = QButtonGroup(self)
        self.tresp_radio = QRadioButton("TrEsp")
        self.tresp_radio.setChecked(True)
        self.tresp_radio.toggled.connect(self.on_parameters_changed)
        self.method_group.addButton(self.tresp_radio)
        settings_layout.addWidget(self.tresp_radio)

        self.dipole_radio = QRadioButton("Dipole")
        self.dipole_radio.toggled.connect(self.on_parameters_changed)
        self.method_group.addButton(self.dipole_radio)
        settings_layout.addWidget(self.dipole_radio)

        cutoff_layout = QHBoxLayout()
        cutoff_layout.addWidget(QLabel("CDC Cutoff"))
        self.cdc_cutoff = QDoubleSpinBox()
        self.cdc_cutoff.setRange(0.0, 1000.0)
        self.cdc_cutoff.setValue(100.0)
        self.cdc_cutoff.setDecimals(1)
        self.cdc_cutoff.setSuffix(" cm⁻¹")
        self.cdc_cutoff.valueChanged.connect(self.on_parameters_changed)
        cutoff_layout.addWidget(self.cdc_cutoff)
        settings_layout.addLayout(cutoff_layout)

        layout.addWidget(settings_group, stretch=0)

        # Site Energies Header
        header_layout = QHBoxLayout()
        site_label = QLabel("Site Energies (Domain cutoff: 20 cm⁻¹)")
        site_label.setProperty("class", "card-title")
        header_layout.addWidget(site_label)
        header_layout.addStretch()

        self.run_hamiltonian_btn = QPushButton("▶ Run Hamiltonian")
        self.run_hamiltonian_btn.clicked.connect(self.run_hamiltonian_requested.emit)
        self.run_hamiltonian_btn.setEnabled(False)
        self.run_hamiltonian_btn.setProperty("class", "primary")
        header_layout.addWidget(self.run_hamiltonian_btn)
        layout.addLayout(header_layout, stretch=0)

        self.site_energy_table = QTableWidget()
        self.site_energy_table.setColumnCount(4)
        self.site_energy_table.setHorizontalHeaderLabels([
            "Pigment", "Residue", "Site Energy (cm⁻¹)", "Domain"
        ])
        header_se = self.site_energy_table.horizontalHeader()
        header_se.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        header_se.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        header_se.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        header_se.setSectionResizeMode(3, QHeaderView.Stretch)
        self.site_energy_table.setAlternatingRowColors(True)
        self.site_energy_table.setEditTriggers(
            QTableWidget.DoubleClicked | QTableWidget.SelectedClicked | QTableWidget.EditKeyPressed
        )
        # Removed inline stylesheet to use global styles
        self.site_energy_table.itemChanged.connect(self.on_site_energy_item_changed)
        layout.addWidget(self.site_energy_table, stretch=1)

        return widget

    def get_parameters(self) -> Dict[str, Any]:
        """Return Hamiltonian-specific parameters"""
        return {
            'coupling_method': 'tresp' if self.tresp_radio.isChecked() else 'dipole',
            'cdc_cutoff': self.cdc_cutoff.value(),
        }

    def on_parameters_changed(self):
        """Emit parameter changes for Hamiltonian settings"""
        self.parameters_changed.emit(self.get_parameters())

    def set_run_enabled(self, enabled: bool):
        """Enable or disable the Run Hamiltonian button"""
        self.run_hamiltonian_btn.setEnabled(enabled)

    def set_hamiltonian_calculator(self, calculator):
        """Set the hamiltonian calculator reference"""
        self.hamiltonian_calculator = calculator

    def update_hamiltonian(self, hamiltonian: np.ndarray, pigment_labels: List[str]):
        """
        Update with Hamiltonian matrix data

        Args:
            hamiltonian: NxN Hamiltonian matrix
            pigment_labels: List of pigment labels
        """
        # Ensure hamiltonian is a numpy array (not DataFrame)
        if hasattr(hamiltonian, 'values'):
            self.hamiltonian = hamiltonian.values
        else:
            self.hamiltonian = np.asarray(hamiltonian)

        self.pigment_labels = pigment_labels

        # Update site energy shift plot
        self.plot_site_energy_shifts()

        # Update info label
        min_coupling = np.min(self.hamiltonian[~np.eye(self.hamiltonian.shape[0], dtype=bool)])
        max_coupling = np.max(self.hamiltonian[~np.eye(self.hamiltonian.shape[0], dtype=bool)])

        self.plot_info_label.setText(
            f"Size: {self.hamiltonian.shape[0]}×{self.hamiltonian.shape[0]} | "
            f"Coupling Range: {min_coupling:.1f} to {max_coupling:.1f} cm⁻¹"
        )
        self.plot_info_label.setStyleSheet("color: #111111;")

        # Enable buttons
        self.view_matrix_btn.setEnabled(True)
        self.export_matrix_btn.setEnabled(True)

        # Update domains using the default cutoff (20 cm⁻¹)
        self.update_domains()

    def update_eigenvalues(self, eigenvalues: np.ndarray, eigenvectors: np.ndarray):
        """
        Update with eigenvalue data

        Args:
            eigenvalues: Array of eigenvalues
            eigenvectors: Matrix of eigenvectors
        """
        self.eigenvalues = eigenvalues
        self.eigenvectors = eigenvectors

    def plot_site_energy_shifts(self):
        """Plot site energy shifts as bar chart"""
        if self.hamiltonian is None:
            return

        self.site_energy_figure.clear()
        ax = self.site_energy_figure.add_subplot(111)

        # Get site energies from hamiltonian diagonal
        site_energies_values = np.diag(self.hamiltonian)

        # Calculate shifts relative to E0a
        shifts = site_energies_values - self.E0a

        # Create bar colors based on shift direction
        colors = ['#d62728' if s > 0 else '#1f77b4' for s in shifts]

        # Create x positions and labels
        x_pos = np.arange(len(self.pigment_labels))
        # Simplify labels (extract just the residue number)
        labels = [name.split('_')[-1] for name in self.pigment_labels]

        bars = ax.bar(x_pos, shifts, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)

        # Add horizontal line at zero
        ax.axhline(y=0, color='black', linestyle='-', linewidth=1)

        # Add value labels on bars (only if not too many bars)
        if len(shifts) <= 25:
            for i, (bar, shift) in enumerate(zip(bars, shifts)):
                height = bar.get_height()
                va = 'bottom' if height >= 0 else 'top'
                offset = 3 if height >= 0 else -3
                ax.annotate(f'{shift:.0f}',
                            xy=(bar.get_x() + bar.get_width() / 2, height),
                            xytext=(0, offset),
                            textcoords="offset points",
                            ha='center', va=va, fontsize=8)

        ax.set_xlabel('Pigment (residue number)', fontsize=10, color='black')
        ax.set_ylabel(f'Site Energy Shift (cm⁻¹) relative to E₀ = {self.E0a} cm⁻¹', fontsize=10, color='black')
        ax.set_title('Site Energy Shifts from CDC Calculation', fontsize=12, color='black', fontweight='bold')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
        ax.grid(axis='y', alpha=0.3)
        ax.tick_params(colors='black')

        # Add statistics annotation
        stats_text = f'Mean shift: {np.mean(shifts):.1f} cm⁻¹\nStd: {np.std(shifts):.1f} cm⁻¹\nRange: {np.min(shifts):.0f} to {np.max(shifts):.0f} cm⁻¹'
        ax.text(0.98, 0.98, stats_text, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        self.site_energy_figure.tight_layout()
        self.site_energy_canvas.draw()


    def update_site_energies_table(self, site_energies: Dict[str, float], pigment_labels: Dict[str, str]):
        """
        Update site energies table

        Args:
            site_energies: Dictionary of {pigment_id: energy}
            pigment_labels: Dictionary of {pigment_id: residue_name}
        """
        # Store site energies for plotting
        self.site_energies = dict(site_energies)

        self._updating_site_table = True
        self.site_energy_table.blockSignals(True)
        self.site_energy_table.setRowCount(len(site_energies))

        for i, (pig_id, energy) in enumerate(site_energies.items()):
            pig_item = QTableWidgetItem(pig_id)
            res_item = QTableWidgetItem(pigment_labels.get(pig_id, ""))
            energy_item = QTableWidgetItem(f"{energy:.1f}")
            energy_item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            energy_item.setToolTip("Double-click to edit site energy")
            energy_item.setData(Qt.UserRole, energy)

            pig_item.setFlags(pig_item.flags() & ~Qt.ItemIsEditable)
            res_item.setFlags(res_item.flags() & ~Qt.ItemIsEditable)
            energy_item.setFlags(energy_item.flags() | Qt.ItemIsEditable)

            self.site_energy_table.setItem(i, 0, pig_item)
            self.site_energy_table.setItem(i, 1, res_item)
            self.site_energy_table.setItem(i, 2, energy_item)
            domain_item = QTableWidgetItem("N/A")
            domain_item.setTextAlignment(Qt.AlignCenter)
            domain_item.setFlags(domain_item.flags() & ~Qt.ItemIsEditable)
            self.site_energy_table.setItem(i, 3, domain_item)

        self.site_energy_table.blockSignals(False)
        self._updating_site_table = False

        # Update domains if Hamiltonian is available
        if self.hamiltonian is not None:
            self.update_domains()

    def on_site_energy_item_changed(self, item: QTableWidgetItem):
        """Handle edits to site energy values"""
        if self._updating_site_table or item.column() != 2:
            return

        pig_item = self.site_energy_table.item(item.row(), 0)
        if pig_item is None:
            return

        pig_id = pig_item.text()
        try:
            energy = float(item.text())
        except (TypeError, ValueError):
            energy = self.site_energies.get(pig_id)
            if energy is None:
                return
            self.site_energy_table.blockSignals(True)
            item.setText(f"{energy:.1f}")
            self.site_energy_table.blockSignals(False)
            return

        self.site_energies[pig_id] = energy
        self.site_energy_table.blockSignals(True)
        item.setText(f"{energy:.1f}")
        item.setData(Qt.UserRole, energy)
        self.site_energy_table.blockSignals(False)
        self.site_energy_updated.emit(pig_id, energy)

    def update_domains(self):
        """Update domain clustering based on cutoff"""
        if self.hamiltonian is None:
            return

        # Ensure hamiltonian is numpy array (safety check for legacy data)
        if hasattr(self.hamiltonian, 'values'):
            self.hamiltonian = self.hamiltonian.values

        # Use fixed cutoff of 20 cm⁻¹
        cutoff = 20.0

        try:
            # Use build_domains method from HamiltonianCalculator if available
            if self.hamiltonian_calculator is not None:
                domains_dict = self.hamiltonian_calculator.build_domains(
                    cutoff=cutoff,
                    hamiltonian=self.hamiltonian,
                    use_absolute=True
                )
                # Convert dict format to list format for compatibility
                domains = [domains_dict[i] for i in sorted(domains_dict.keys())]
            else:
                # Fall back to local cluster_domains method
                domains = self.cluster_domains(self.hamiltonian, cutoff)
                # Convert list format to dict format for emission
                domains_dict = {i: domain for i, domain in enumerate(domains)}
        except Exception as e:
            print(f"Error in domain clustering: {e}")
            print(f"Hamiltonian type: {type(self.hamiltonian)}")
            print(f"Hamiltonian shape: {self.hamiltonian.shape if hasattr(self.hamiltonian, 'shape') else 'no shape'}")
            import traceback
            traceback.print_exc()
            raise

        # Update site energy table domain column
        self.update_site_energy_domains(domains)

        # Emit signal with domains dictionary for 3D visualization
        self.domains_updated.emit(domains_dict)

    def cluster_domains(self, hamiltonian, cutoff):
        """
        Find connected components in coupling graph

        Args:
            hamiltonian: Hamiltonian matrix
            cutoff: Coupling cutoff threshold (cm⁻¹)

        Returns:
            List of domains, where each domain is a list of pigment indices
        """
        # Ensure hamiltonian is a numpy array
        if hasattr(hamiltonian, 'values'):
            hamiltonian = hamiltonian.values
        else:
            hamiltonian = np.asarray(hamiltonian)

        n = hamiltonian.shape[0]
        graph = {i: [] for i in range(n)}

        # Build adjacency graph (off-diagonal couplings above cutoff)
        for i in range(n):
            for j in range(i+1, n):
                if abs(hamiltonian[i, j]) >= cutoff:
                    graph[i].append(j)
                    graph[j].append(i)

        # Find connected components using DFS
        visited = set()
        domains = []
        for node in range(n):
            if node not in visited:
                domain = []
                stack = [node]
                while stack:
                    current = stack.pop()
                    if current not in visited:
                        visited.add(current)
                        domain.append(current)
                        stack.extend(graph[current])
                domains.append(sorted(domain))

        return domains

    def update_site_energy_domains(self, domains):
        """Update domain assignments in site energy table"""
        if not self.pigment_labels:
            return

        # Create mapping of pigment index to domain
        domain_map = {}
        for domain_idx, domain in enumerate(domains):
            for pig_idx in domain:
                domain_map[pig_idx] = domain_idx + 1

        # Update table - match pigment IDs to get correct indices
        for row in range(self.site_energy_table.rowCount()):
            # Get pigment ID from table
            pig_id_item = self.site_energy_table.item(row, 0)
            if pig_id_item:
                pig_id = pig_id_item.text()
                # Find index in pigment_labels
                try:
                    pig_idx = self.pigment_labels.index(pig_id)
                    domain_id = domain_map.get(pig_idx, 0)
                except ValueError:
                    domain_id = 0
            else:
                domain_id = 0

            domain_item = QTableWidgetItem(f"D{domain_id}" if domain_id > 0 else "N/A")
            domain_item.setTextAlignment(Qt.AlignCenter)
            domain_item.setFlags(domain_item.flags() & ~Qt.ItemIsEditable)
            self.site_energy_table.setItem(row, 3, domain_item)

    def get_current_domains(self):
        """
        Get current domain clustering for spectra calculation

        Returns:
            dict: {domain_id: [pigment_indices]}
        """
        if self.hamiltonian is None:
            return {}

        # Use fixed cutoff of 20 cm⁻¹
        cutoff = 20.0

        # Use build_domains method from HamiltonianCalculator if available
        if self.hamiltonian_calculator is not None:
            domains_dict = self.hamiltonian_calculator.build_domains(
                cutoff=cutoff,
                hamiltonian=self.hamiltonian,
                use_absolute=True
            )
        else:
            # Fall back to local cluster_domains method
            domains_list = self.cluster_domains(self.hamiltonian, cutoff)
            # Convert to dict format
            domains_dict = {}
            for domain_id, pigment_indices in enumerate(domains_list):
                domains_dict[domain_id] = pigment_indices

        return domains_dict

    def clear(self):
        """Clear all data"""
        self.hamiltonian = None
        self.eigenvalues = None
        self.eigenvectors = None
        self.pigment_labels = None
        self.site_energies = None

        self.site_energy_figure.clear()
        self.site_energy_canvas.draw()
        self.site_energy_table.setRowCount(0)

        self.plot_info_label.setText("Site energies not calculated")
        self.plot_info_label.setStyleSheet("color: #6b6b6b; font-style: italic;")

        self.view_matrix_btn.setEnabled(False)
        self.export_matrix_btn.setEnabled(False)
