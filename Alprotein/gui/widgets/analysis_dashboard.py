"""
Analysis Dashboard Widget for Scientific Workbench

Provides a three-panel visualization dashboard showing:
- Site Energies plot
- Hamiltonian heatmap
- Absorption Spectrum
"""

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QFrame, QSizePolicy, QTabWidget
)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
from typing import Dict, List, Optional


class MiniPlotWidget(QWidget):
    """Base class for mini plot widgets in the dashboard"""

    export_requested = pyqtSignal()
    view_details = pyqtSignal()

    def __init__(self, title, parent=None):
        super().__init__(parent)
        self.title = title
        self.setup_ui()

    def setup_ui(self):
        """Setup the UI"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(6)

        # Header
        header_layout = QHBoxLayout()

        title_label = QLabel(self.title)
        title_label.setStyleSheet("""
            QLabel {
                font-weight: bold;
                color: #111111;
                font-size: 11px;
            }
        """)
        header_layout.addWidget(title_label)

        header_layout.addStretch()

        # Action buttons
        self.details_btn = QPushButton("Details")
        self.details_btn.setMaximumWidth(70)
        self.details_btn.setStyleSheet("""
            QPushButton {
                background-color: #111111;
                color: #ffffff;
                border: 1px solid #111111;
                padding: 4px 8px;
                border-radius: 2px;
                font-size: 10px;
            }
            QPushButton:hover {
                background-color: #000000;
            }
        """)
        self.details_btn.clicked.connect(self.view_details.emit)
        header_layout.addWidget(self.details_btn)

        self.export_btn = QPushButton("Export")
        self.export_btn.setMaximumWidth(70)
        self.export_btn.setStyleSheet("""
            QPushButton {
                background-color: #ffffff;
                color: #111111;
                border: 1px solid #111111;
                padding: 4px 8px;
                border-radius: 2px;
                font-size: 10px;
            }
            QPushButton:hover {
                background-color: #efefef;
            }
        """)
        self.export_btn.clicked.connect(self.export_requested.emit)
        header_layout.addWidget(self.export_btn)

        layout.addLayout(header_layout)

        # Plot area - Enhanced with larger size and higher DPI
        self.figure = Figure(figsize=(7, 5), facecolor='white', dpi=100)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.canvas.setMinimumSize(500, 400)  # Set minimum size for better visibility
        layout.addWidget(self.canvas)

    def clear_plot(self):
        """Clear the plot"""
        self.figure.clear()
        self.canvas.draw()

    def plot_histogram(self, data: np.ndarray, xlabel: str = "", title: str = ""):
        """Plot a histogram"""
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.set_facecolor('white')

        ax.hist(data, bins=30, color='#2b2b2b', alpha=0.7, edgecolor='black', linewidth=0.8)

        if xlabel:
            ax.set_xlabel(xlabel, fontsize=9, color='black')
        ax.set_ylabel('Count', fontsize=9, color='black')
        if title:
            ax.set_title(title, fontsize=10, color='black', fontweight='bold')

        ax.grid(axis='y', alpha=0.3, linestyle='--')
        ax.tick_params(axis='both', labelsize=8, colors='black')

        self.figure.tight_layout()
        self.canvas.draw()

    def plot_heatmap(self, matrix: np.ndarray, title: str = ""):
        """Plot a heatmap"""
        self.figure.clear()
        ax = self.figure.add_subplot(111)

        # Plot heatmap
        im = ax.imshow(matrix, cmap='RdBu_r', aspect='auto', interpolation='nearest')

        if title:
            ax.set_title(title, fontsize=10, color='black', fontweight='bold')

        # Minimal labeling for space
        n = matrix.shape[0]
        if n <= 10:
            ax.set_xticks(range(n))
            ax.set_yticks(range(n))
            ax.tick_params(axis='both', labelsize=7, colors='black')
        else:
            ax.set_xticks([])
            ax.set_yticks([])

        self.figure.tight_layout()
        self.canvas.draw()


class SiteEnergyMiniPlot(MiniPlotWidget):
    """Mini site energy plot"""

    def __init__(self, parent=None):
        super().__init__("Site Energies", parent)

    def update_plot(self, site_energies: Dict[str, float], vacuum_energies: Dict[str, float]):
        """Update with site energy data"""
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.set_facecolor('white')

        pigment_ids = list(site_energies.keys())
        energies = [site_energies[pid] for pid in pigment_ids]
        vacuums = [vacuum_energies.get(pid, 0) for pid in pigment_ids]
        shifts = [e - v for e, v in zip(energies, vacuums)]

        x_pos = np.arange(len(pigment_ids))

        # Plot as scatter with color coding
        scatter = ax.scatter(x_pos, energies, s=60, c=shifts,
                           cmap='RdYlGn_r', alpha=0.8, edgecolors='black',
                           linewidth=1, zorder=3)

        # Add stems
        for i, (x, y) in enumerate(zip(x_pos, energies)):
            ax.plot([x, x], [vacuums[i], y], 'k-', alpha=0.2, linewidth=1)

        ax.set_ylabel('Energy (cm⁻¹)', fontsize=9, color='black')
        ax.set_xticks([])
        ax.grid(axis='y', alpha=0.2)
        ax.tick_params(axis='both', labelsize=8, colors='black')

        self.figure.tight_layout()
        self.canvas.draw()


class HamiltonianMiniPlot(MiniPlotWidget):
    """Mini Hamiltonian heatmap"""

    def __init__(self, parent=None):
        super().__init__("Hamiltonian", parent)

    def update_plot(self, hamiltonian: np.ndarray):
        """Update with Hamiltonian matrix"""
        self.figure.clear()
        ax = self.figure.add_subplot(111)

        im = ax.imshow(hamiltonian, cmap='RdBu_r', aspect='auto',
                      interpolation='nearest')

        # Minimal labeling for space
        n = hamiltonian.shape[0]
        if n <= 10:
            ax.set_xticks(range(n))
            ax.set_yticks(range(n))
            ax.tick_params(axis='both', labelsize=7, colors='black')
        else:
            ax.set_xticks([])
            ax.set_yticks([])

        self.figure.tight_layout()
        self.canvas.draw()


class SpectrumMiniPlot(MiniPlotWidget):
    """Mini absorption spectrum plot"""

    def __init__(self, parent=None):
        super().__init__("Spectrum", parent)

    def update_plot(self, wavelengths: np.ndarray, spectrum: np.ndarray):
        """Update with spectrum data"""
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.set_facecolor('white')

        ax.plot(wavelengths, spectrum, 'b-', linewidth=2)
        ax.fill_between(wavelengths, spectrum, alpha=0.3, color='blue')

        ax.set_xlabel('λ (nm)', fontsize=9, color='black')
        ax.set_ylabel('Intensity', fontsize=9, color='black')
        ax.grid(True, alpha=0.2)
        ax.tick_params(axis='both', labelsize=8, colors='black')

        self.figure.tight_layout()
        self.canvas.draw()


class AnalysisDashboard(QFrame):
    """
    Analysis dashboard showing three mini plots side-by-side
    """

    view_site_energies = pyqtSignal()
    view_hamiltonian = pyqtSignal()
    view_spectrum = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setup_ui()

    def setup_ui(self):
        """Setup the UI"""
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(16, 16, 16, 16)
        main_layout.setSpacing(16)

        # Header
        header = QLabel("ANALYSIS DASHBOARD")
        header.setStyleSheet("""
            QLabel {
                font-size: 13px;
                font-weight: bold;
                color: #111111;
                padding: 8px;
                background-color: #f3f3f3;
                border-radius: 2px;
            }
        """)
        main_layout.addWidget(header)

        # Three plots in horizontal layout
        plots_layout = QHBoxLayout()
        plots_layout.setSpacing(16)

        # Site Energy plot
        self.site_energy_plot = SiteEnergyMiniPlot()
        self.site_energy_plot.view_details.connect(self.view_site_energies.emit)
        plots_layout.addWidget(self.site_energy_plot, stretch=1)

        # Hamiltonian plot
        self.hamiltonian_plot = HamiltonianMiniPlot()
        self.hamiltonian_plot.view_details.connect(self.view_hamiltonian.emit)
        plots_layout.addWidget(self.hamiltonian_plot, stretch=1)

        # Spectrum plot
        self.spectrum_plot = SpectrumMiniPlot()
        self.spectrum_plot.view_details.connect(self.view_spectrum.emit)
        plots_layout.addWidget(self.spectrum_plot, stretch=1)

        main_layout.addLayout(plots_layout)

        # Set frame style
        self.setFrameShape(QFrame.StyledPanel)
        self.setStyleSheet("""
            AnalysisDashboard {
                background-color: #ffffff;
                border: 1px solid #e1e1e1;
                border-radius: 2px;
            }
        """)

    def update_site_energies(self, site_energies: Dict[str, float], vacuum_energies: Dict[str, float]):
        """Update site energy plot"""
        self.site_energy_plot.update_plot(site_energies, vacuum_energies)

    def update_hamiltonian(self, hamiltonian: np.ndarray):
        """Update Hamiltonian plot"""
        self.hamiltonian_plot.update_plot(hamiltonian)

    def update_spectrum(self, wavelengths: np.ndarray, spectrum: np.ndarray):
        """Update spectrum plot"""
        self.spectrum_plot.update_plot(wavelengths, spectrum)

    def clear_all(self):
        """Clear all plots"""
        self.site_energy_plot.clear_plot()
        self.hamiltonian_plot.clear_plot()
        self.spectrum_plot.clear_plot()
