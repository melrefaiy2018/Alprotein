"""
Enhanced Results Panel Widget for Alprotein GUI

This widget provides a comprehensive multi-tab results display:
- Overview: Summary statistics
- Site Energy: Energy level diagrams with hover details
- Hamiltonian: Matrix visualization and eigenvalues
- Eigenvalues: Detailed table
- Absorption Spectrum: Publication-quality plot
- CDC Analysis: Contribution breakdown
"""

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTabWidget, QLabel,
    QPushButton, QTableWidget, QTableWidgetItem, QHeaderView,
    QGroupBox, QSplitter, QTextEdit, QSizePolicy
)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont, QColor, QBrush

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
from typing import Dict, Any, List, Optional

# Import our custom widgets
from .spectrum_widget import SpectrumPlotWidget
from .hamiltonian_widget import HamiltonianWidget
from .results_panel import ContributionAnalysisWidget


class SiteEnergyPlotWidget(QWidget):
    """Widget for plotting site energies as a scatter/bar plot"""

    pigment_clicked = pyqtSignal(str)  # pigment_id

    def __init__(self, parent=None):
        super().__init__(parent)
        self.site_energies = {}
        self.vacuum_energies = {}
        self.setup_ui()

    def setup_ui(self):
        """Setup the UI"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)

        # Create matplotlib figure
        self.figure = Figure(figsize=(10, 5), facecolor='white', dpi=100)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        layout.addWidget(self.canvas)

        # Info label
        self.info_label = QLabel("Hover over points for details")
        self.info_label.setStyleSheet("color: #6b6b6b; font-style: italic;")
        self.info_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(self.info_label)

        # Export buttons
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()

        self.export_csv_btn = QPushButton("📊 Export CSV")
        self.export_csv_btn.clicked.connect(self.export_csv)
        self.export_csv_btn.setStyleSheet("""
            QPushButton {
                background-color: #111111;
                color: #ffffff;
                border: 1px solid #111111;
                padding: 6px 12px;
                border-radius: 2px;
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
                border: 1px solid #111111;
                padding: 6px 12px;
                border-radius: 2px;
            }
            QPushButton:hover {
                background-color: #efefef;
            }
        """)
        btn_layout.addWidget(self.export_plot_btn)

        layout.addLayout(btn_layout)

    def update_data(self, site_energies: Dict[str, float], vacuum_energies: Dict[str, float]):
        """Update with site energy data"""
        self.site_energies = site_energies
        self.vacuum_energies = vacuum_energies
        self.plot_energies()

    def plot_energies(self):
        """Plot site energies"""
        if not self.site_energies:
            return

        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.set_facecolor('white')

        # Prepare data
        pigment_ids = list(self.site_energies.keys())
        site_energies = [self.site_energies[pid] for pid in pigment_ids]
        vacuum_energies = [self.vacuum_energies.get(pid, 0) for pid in pigment_ids]
        energy_shifts = [se - ve for se, ve in zip(site_energies, vacuum_energies)]

        x_pos = np.arange(len(pigment_ids))

        # Plot vacuum energies as a baseline
        ax.axhline(y=np.mean(vacuum_energies), color='gray', linestyle='--',
                  alpha=0.5, label='Average Vacuum Energy')

        # Plot site energies as scatter
        scatter = ax.scatter(x_pos, site_energies, s=100, c=energy_shifts,
                           cmap='RdYlGn_r', alpha=0.8, edgecolors='black',
                           linewidth=1.5, zorder=3)

        # Add stems
        for i, (x, y) in enumerate(zip(x_pos, site_energies)):
            ax.plot([x, x], [vacuum_energies[i], y], 'k-', alpha=0.3, linewidth=1)

        # Colorbar
        cbar = self.figure.colorbar(scatter, ax=ax)
        cbar.set_label('Energy Shift (cm⁻¹)', rotation=270, labelpad=20, color='black')
        cbar.ax.tick_params(colors='black')

        # Labels and formatting
        ax.set_xlabel('Pigment ID', fontsize=12, color='black')
        ax.set_ylabel('Site Energy (cm⁻¹)', fontsize=12, color='black')
        ax.set_title('Site Energies of Pigments', fontsize=14, color='black', fontweight='bold')

        ax.set_xticks(x_pos)
        ax.set_xticklabels(pigment_ids, rotation=45, ha='right', fontsize=9)

        ax.grid(axis='y', alpha=0.3, linestyle='--')
        ax.legend(loc='best')

        # Set tick colors
        ax.tick_params(colors='black')
        for spine in ax.spines.values():
            spine.set_edgecolor('black')

        self.figure.tight_layout()
        self.canvas.draw()

        # Update info
        mean_shift = np.mean(energy_shifts)
        std_shift = np.std(energy_shifts)
        self.info_label.setText(
            f"Mean shift: {mean_shift:+.1f} cm⁻¹ | Std dev: {std_shift:.1f} cm⁻¹"
        )
        self.info_label.setStyleSheet("color: #111111; background-color: transparent; font-weight: bold;")

    def export_csv(self):
        """Export site energies to CSV"""
        from PyQt5.QtWidgets import QFileDialog, QMessageBox
        import sys
        import csv

        if not self.site_energies:
            return

        options = QFileDialog.DontUseNativeDialog if sys.platform == 'darwin' else QFileDialog.Options()

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Site Energies",
            "site_energies.csv",
            "CSV Files (*.csv);;All Files (*)",
            options=options
        )

        if file_path:
            try:
                with open(file_path, 'w', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow(['Pigment', 'Site Energy (cm-1)', 'Vacuum Energy (cm-1)', 'Energy Shift (cm-1)'])

                    for pigment_id in self.site_energies.keys():
                        site_energy = self.site_energies[pigment_id]
                        vacuum_energy = self.vacuum_energies.get(pigment_id, 0)
                        energy_shift = site_energy - vacuum_energy
                        writer.writerow([pigment_id, f"{site_energy:.1f}",
                                       f"{vacuum_energy:.1f}", f"{energy_shift:+.1f}"])

                QMessageBox.information(self, "Export Successful",
                                      f"Site energies exported to:\n{file_path}")
            except Exception as e:
                QMessageBox.critical(self, "Export Error",
                                   f"Failed to export:\n{str(e)}")

    def export_plot(self):
        """Export plot as image"""
        from PyQt5.QtWidgets import QFileDialog
        import sys

        options = QFileDialog.DontUseNativeDialog if sys.platform == 'darwin' else QFileDialog.Options()

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Plot",
            "site_energies_plot.png",
            "PNG Files (*.png);;PDF Files (*.pdf);;All Files (*)",
            options=options
        )

        if file_path:
            self.figure.savefig(file_path, dpi=300, bbox_inches='tight',
                              facecolor='white', edgecolor='none')

    def clear(self):
        """Clear the plot"""
        self.site_energies = {}
        self.vacuum_energies = {}
        self.figure.clear()
        self.canvas.draw()
        self.info_label.setText("Hover over points for details")
        self.info_label.setStyleSheet("color: #6b6b6b; font-style: italic;")


class OverviewWidget(QWidget):
    """Overview summary widget"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setup_ui()

    def setup_ui(self):
        """Setup the UI"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(15)

        # Title
        title = QLabel("CALCULATION SUMMARY")
        title.setAlignment(Qt.AlignCenter)
        title.setStyleSheet("""
            QLabel {
                background-color: #f3f3f3;
                color: #111111;
                padding: 15px;
                font-size: 16px;
                font-weight: bold;
                letter-spacing: 2px;
                border-radius: 2px;
                border: 1px solid #e1e1e1;
            }
        """)
        layout.addWidget(title)

        # Summary text
        self.summary_text = QTextEdit()
        self.summary_text.setReadOnly(True)
        self.summary_text.setStyleSheet("""
            QTextEdit {
                background-color: #ffffff;
                color: #111111;
                border: 1px solid #e1e1e1;
                border-radius: 2px;
                padding: 15px;
                font-size: 12px;
                font-family: 'Courier New', monospace;
            }
        """)
        self.summary_text.setPlainText("No results available yet.\n\nLoad a structure and run calculations to see summary.")
        layout.addWidget(self.summary_text)

    def update_summary(self, results: Dict[str, Any]):
        """Update summary with results"""
        summary = "=" * 60 + "\n"
        summary += "ALPROTEIN CALCULATION SUMMARY\n"
        summary += "=" * 60 + "\n\n"

        # Site energies
        if 'site_energies' in results:
            se_data = results['site_energies']
            site_energies = se_data.get('site_energies', {})
            vacuum_energies = se_data.get('vacuum_energies', {})

            summary += "SITE ENERGIES\n"
            summary += "-" * 40 + "\n"
            summary += f"Number of pigments: {len(site_energies)}\n"

            if site_energies:
                energies = list(site_energies.values())
                vacuums = [vacuum_energies.get(pid, 0) for pid in site_energies.keys()]
                shifts = [e - v for e, v in zip(energies, vacuums)]

                summary += f"Energy range: {min(energies):.1f} - {max(energies):.1f} cm⁻¹\n"
                summary += f"Mean site energy: {np.mean(energies):.1f} cm⁻¹\n"
                summary += f"Mean energy shift: {np.mean(shifts):+.1f} cm⁻¹\n"
                summary += f"Std dev shift: {np.std(shifts):.1f} cm⁻¹\n"

            summary += "\n"

        # Hamiltonian
        if 'hamiltonian' in results:
            ham_data = results['hamiltonian']
            hamiltonian = ham_data.get('matrix')
            eigenvalues = ham_data.get('eigenvalues')

            if hamiltonian is not None:
                summary += "HAMILTONIAN\n"
                summary += "-" * 40 + "\n"
                summary += f"Matrix size: {hamiltonian.shape[0]}×{hamiltonian.shape[1]}\n"

                # Get off-diagonal elements (couplings)
                mask = ~np.eye(hamiltonian.shape[0], dtype=bool)
                couplings = hamiltonian[mask]

                summary += f"Coupling range: {np.min(couplings):.1f} - {np.max(couplings):.1f} cm⁻¹\n"
                summary += f"Mean coupling: {np.mean(couplings):.1f} cm⁻¹\n"

                if eigenvalues is not None:
                    summary += f"\nEigenvalue range: {np.min(eigenvalues):.1f} - {np.max(eigenvalues):.1f} cm⁻¹\n"
                    summary += f"Energy spread: {np.max(eigenvalues) - np.min(eigenvalues):.1f} cm⁻¹\n"

                summary += "\n"

        # Spectrum
        if 'spectrum' in results:
            spec_data = results['spectrum']
            summary += "ABSORPTION SPECTRUM\n"
            summary += "-" * 40 + "\n"
            summary += f"Peak wavelength: {spec_data.get('peak_wavelength', 0):.1f} nm\n"
            summary += f"Peak energy: {spec_data.get('peak_energy', 0):.1f} cm⁻¹\n"
            summary += f"FWHM: {spec_data.get('fwhm', 0):.1f} nm\n"
            summary += "\n"

        summary += "=" * 60 + "\n"

        self.summary_text.setPlainText(summary)


class EnhancedResultsPanel(QWidget):
    """
    Enhanced results panel with multi-tab display
    """

    pigment_selected = pyqtSignal(str)
    export_requested = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.results = {}
        self.setup_ui()

    def setup_ui(self):
        """Setup the UI"""
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # Create tab widget
        self.tab_widget = QTabWidget()
        self.tab_widget.setStyleSheet("""
            QTabWidget::pane {
                background-color: #ffffff;
                border: 1px solid #e1e1e1;
            }
            QTabBar::tab {
                color: #111111;
                background-color: #f3f3f3;
                padding: 10px 18px;
                margin: 2px;
                border: 1px solid #e1e1e1;
                border-bottom: none;
                border-top-left-radius: 2px;
                border-top-right-radius: 2px;
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

        # Overview tab
        self.overview_widget = OverviewWidget()
        self.tab_widget.addTab(self.overview_widget, "📊 Overview")

        # Site Energy tab
        self.site_energy_widget = SiteEnergyPlotWidget()
        self.tab_widget.addTab(self.site_energy_widget, "⚡ Site Energy")

        # Hamiltonian tab
        self.hamiltonian_widget = HamiltonianWidget()
        self.tab_widget.addTab(self.hamiltonian_widget, "🔢 Hamiltonian")

        # Absorption Spectrum tab
        self.spectrum_widget = SpectrumPlotWidget()
        self.tab_widget.addTab(self.spectrum_widget, "📈 Absorption Spectrum")

        # CDC Analysis tab
        self.cdc_widget = ContributionAnalysisWidget()
        self.cdc_widget.pigment_selected.connect(self.pigment_selected.emit)
        self.tab_widget.addTab(self.cdc_widget, "🔬 CDC Analysis")

        main_layout.addWidget(self.tab_widget)

    def display_site_energies(self, results: Dict[str, Any]):
        """Display site energy results"""
        site_energies = results.get('site_energies', {})
        vacuum_energies = results.get('vacuum_energies', {})

        self.site_energy_widget.update_data(site_energies, vacuum_energies)

        # Update overview
        self.results['site_energies'] = results
        self.overview_widget.update_summary(self.results)

    def display_hamiltonian(self, hamiltonian: np.ndarray, pigment_labels: List[str]):
        """Display Hamiltonian matrix"""
        self.hamiltonian_widget.update_hamiltonian(hamiltonian, pigment_labels)

        # Store in results
        self.results['hamiltonian'] = {
            'matrix': hamiltonian,
            'pigment_labels': pigment_labels
        }
        self.overview_widget.update_summary(self.results)

    def display_eigenvalues(self, eigenvalues: np.ndarray, eigenvectors: np.ndarray):
        """Display eigenvalues and eigenvectors"""
        self.hamiltonian_widget.update_eigenvalues(eigenvalues, eigenvectors)

        # Update results
        if 'hamiltonian' in self.results:
            self.results['hamiltonian']['eigenvalues'] = eigenvalues
            self.results['hamiltonian']['eigenvectors'] = eigenvectors
        else:
            self.results['hamiltonian'] = {
                'eigenvalues': eigenvalues,
                'eigenvectors': eigenvectors
            }

        self.overview_widget.update_summary(self.results)

    def display_spectrum(self, eigenvalues: np.ndarray, oscillator_strengths: np.ndarray):
        """Display absorption spectrum (legacy method for backward compatibility)"""
        self.spectrum_widget.update_spectrum_data(eigenvalues, oscillator_strengths)

        # Store in results
        self.results['spectrum'] = {
            'eigenvalues': eigenvalues,
            'oscillator_strengths': oscillator_strengths
        }
        self.overview_widget.update_summary(self.results)

    def display_spectra(self, wavelengths_abs: np.ndarray, absorption: np.ndarray,
                       wavelengths_fl: np.ndarray, fluorescence: np.ndarray):
        """
        Display both absorption and fluorescence spectra

        Args:
            wavelengths_abs: Wavelength axis for absorption (nm)
            absorption: Absorption spectrum (normalized)
            wavelengths_fl: Wavelength axis for fluorescence (nm)
            fluorescence: Fluorescence spectrum (normalized)
        """
        # Update the spectrum widget with both spectra
        self.spectrum_widget.update_spectra(
            wavelengths_abs, absorption,
            wavelengths_fl, fluorescence
        )

        # Find peaks for summary
        peak_idx_abs = np.argmax(absorption)
        peak_wl_abs = wavelengths_abs[peak_idx_abs]

        peak_idx_fl = np.argmax(fluorescence)
        peak_wl_fl = wavelengths_fl[peak_idx_fl]

        # Store in results
        self.results['spectrum'] = {
            'wavelengths_abs': wavelengths_abs,
            'absorption': absorption,
            'wavelengths_fl': wavelengths_fl,
            'fluorescence': fluorescence,
            'peak_wavelength': peak_wl_abs,
            'peak_energy': 1e7 / peak_wl_abs,
            'fl_peak_wavelength': peak_wl_fl,
            'fl_peak_energy': 1e7 / peak_wl_fl,
            'stokes_shift': 1e7 / peak_wl_abs - 1e7 / peak_wl_fl,
            'fwhm': 0.0  # Will be calculated if needed
        }
        self.overview_widget.update_summary(self.results)

    def display_cdc_analysis(self, results: Dict[str, Any]):
        """Display CDC analysis results"""
        self.cdc_widget.update_data(results.get('detailed_contributions', {}))

        # Store in results
        self.results['cdc_analysis'] = results
        self.overview_widget.update_summary(self.results)

    def clear_results(self):
        """Clear all results"""
        self.results = {}
        self.site_energy_widget.clear()
        self.hamiltonian_widget.clear()
        self.spectrum_widget.clear()
        self.overview_widget.summary_text.setPlainText(
            "No results available yet.\n\nLoad a structure and run calculations to see summary."
        )
