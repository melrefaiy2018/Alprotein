#!/usr/bin/env python3
"""
Test script to verify exciton distribution plotting in GUI.
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')  # Use Qt5 backend

# Add Alprotein to path
sys.path.insert(0, str(Path(__file__).parent))

from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget
from Alprotein.gui.widgets.analysis_dashboard import MiniPlotWidget
import matplotlib.pyplot as plt


def create_test_distribution():
    """Create test exciton distribution data."""
    distributions = {}
    labels = ['A_CLA_501', 'A_CLA_502', 'A_CLA_503']

    for i, label in enumerate(labels):
        # Create realistic distribution (Gaussian-like)
        energies = np.linspace(14500, 15500, 50)
        center = 15000 + i * 100
        probs = np.exp(-((energies - center) ** 2) / (2 * 50 ** 2))
        probs = probs / np.sum(probs)
        distributions[label] = (energies.tolist(), probs.tolist())

    return distributions, labels


def test_miniplot_widget():
    """Test MiniPlotWidget with exciton distribution."""
    app = QApplication(sys.argv)

    # Create main window
    window = QMainWindow()
    window.setWindowTitle("Exciton Distribution Plot Test")
    window.resize(800, 600)

    # Create central widget
    central = QWidget()
    layout = QVBoxLayout(central)

    # Create MiniPlotWidget
    plot_widget = MiniPlotWidget("Test Exciton Distribution")
    layout.addWidget(plot_widget)

    window.setCentralWidget(central)

    # Get test data
    distributions, labels = create_test_distribution()

    # Plot the data (similar to update_exciton_plots)
    fig = plot_widget.figure
    fig.clear()
    ax = fig.add_subplot(111)

    # Plot each pigment's distribution
    colors = plt.cm.tab20(np.linspace(0, 1, len(labels)))
    for i, label in enumerate(labels):
        energies, probs = distributions[label]
        wavelengths = [1e7 / e if e > 0 else 0 for e in energies]
        ax.hist(wavelengths, bins=30, weights=probs, alpha=0.5,
                label=label.split('_')[-1], color=colors[i])

    ax.set_xlabel('Wavelength (nm)', fontsize=10)
    ax.set_ylabel('Probability', fontsize=10)
    ax.set_title('Test Exciton Distribution', fontsize=11, fontweight='bold')
    ax.set_xlim(600, 720)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()

    # Draw the canvas
    plot_widget.canvas.draw()
    plot_widget.canvas.flush_events()
    plot_widget.update()

    print("✓ Plot widget created and data plotted")
    print(f"✓ Figure: {fig}")
    print(f"✓ Canvas: {plot_widget.canvas}")
    print(f"✓ Distribution labels: {labels}")

    # Show the window
    window.show()

    print("\n==> If you see a plot with colored histograms, the plotting works!")
    print("==> Close the window to exit.")

    sys.exit(app.exec_())


if __name__ == "__main__":
    test_miniplot_widget()
