"""
Summary Cards Widget for Data Analysis Tab

Displays three summary cards with key statistics.
"""

from PyQt5.QtWidgets import (
    QWidget, QHBoxLayout, QVBoxLayout, QLabel, QGroupBox
)
from PyQt5.QtCore import Qt
from typing import Dict
import numpy as np


class SummaryCard(QGroupBox):
    """Individual summary card"""

    def __init__(self, title, color="#111111"):
        super().__init__(title)
        self.color = color
        self.setup_ui()

    def setup_ui(self):
        """Setup card UI"""
        self.setProperty("class", "card")
        
        layout = QVBoxLayout(self)
        layout.setSpacing(8)
        layout.setContentsMargins(16, 24, 16, 16)

        self.content_labels = []

    def add_stat(self, name, value):
        """Add a statistic to the card"""
        stat_layout = QHBoxLayout()

        name_label = QLabel(f"{name}:")
        name_label.setProperty("class", "text-secondary")
        stat_layout.addWidget(name_label)

        value_label = QLabel(str(value))
        value_label.setStyleSheet(f"color: {self.color}; font-weight: 600; font-size: 14px;")
        stat_layout.addWidget(value_label)

        stat_layout.addStretch()
        self.layout().addLayout(stat_layout)
        self.content_labels.append((name_label, value_label))

    def clear_stats(self):
        """Clear all statistics"""
        for name_label, value_label in self.content_labels:
            name_label.deleteLater()
            value_label.deleteLater()
        self.content_labels = []


class SummaryCardsWidget(QWidget):
    """
    Widget displaying three summary cards
    """

    def __init__(self):
        super().__init__()
        self.setup_ui()

    def setup_ui(self):
        """Setup the UI"""
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(16)

        # Card 1: Statistical Summary
        self.stats_card = SummaryCard("📊 Statistical Summary", "#111111")
        layout.addWidget(self.stats_card)

        # Card 2: Pigment Counts
        self.pigments_card = SummaryCard("🧬 Pigment Counts", "#333333")
        layout.addWidget(self.pigments_card)

        # Card 3: Interaction Summary
        self.interaction_card = SummaryCard("⚡ Interaction Summary", "#555555")
        layout.addWidget(self.interaction_card)

        self.setMaximumHeight(150)

    def update_statistics(self, site_energies: Dict[str, float]):
        """Update statistical summary card"""
        self.stats_card.clear_stats()

        if site_energies:
            energies = list(site_energies.values())
            self.stats_card.add_stat("Mean", f"{np.mean(energies):.1f} cm⁻¹")
            self.stats_card.add_stat("Std Dev", f"{np.std(energies):.1f} cm⁻¹")
            self.stats_card.add_stat("Min", f"{np.min(energies):.1f} cm⁻¹")
            self.stats_card.add_stat("Max", f"{np.max(energies):.1f} cm⁻¹")

    def update_pigments(self, pigment_counts: Dict[str, int]):
        """Update pigment counts card"""
        self.pigments_card.clear_stats()

        total = sum(pigment_counts.values())
        self.pigments_card.add_stat("Total", str(total))

        for pig_type, count in pigment_counts.items():
            self.pigments_card.add_stat(pig_type, str(count))

    def update_interactions(self, hamiltonian: np.ndarray):
        """Update interaction summary card"""
        self.interaction_card.clear_stats()

        if hamiltonian is not None:
            # Convert to numpy array if it's a DataFrame
            ham_array = hamiltonian.values if hasattr(hamiltonian, 'values') else hamiltonian
            n = ham_array.shape[0]
            off_diag = []
            for i in range(n):
                for j in range(i+1, n):
                    off_diag.append(abs(ham_array[i, j]))

            if off_diag:
                self.interaction_card.add_stat("Min Coupling", f"{min(off_diag):.1f} cm⁻¹")
                self.interaction_card.add_stat("Avg Coupling", f"{np.mean(off_diag):.1f} cm⁻¹")
                self.interaction_card.add_stat("Max Coupling", f"{max(off_diag):.1f} cm⁻¹")
                strong = sum(1 for c in off_diag if c > 100)
                self.interaction_card.add_stat("Strong (>100)", str(strong))
