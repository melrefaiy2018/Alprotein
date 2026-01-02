"""
Actions Panel Widget
"""

from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QPushButton, QLabel, QGroupBox)
from PyQt5.QtCore import Qt, pyqtSignal


class ActionsPanel(QWidget):
    """
    Panel containing action buttons for calculations
    """
    calculate_site_energies = pyqtSignal()
    run_cdc_analysis = pyqtSignal()

    def __init__(self):
        super().__init__()
        self.site_energies_calculated = False
        self.system_loaded = False
        self.setup_ui()
        self.connect_signals()
        self.update_button_states()

    def setup_ui(self):
        """Setup the UI components"""
        layout = QVBoxLayout(self)
        
        # Actions group
        group_box = QGroupBox("🚀 Calculations")
        group_box.setProperty("class", "card")
        group_layout = QVBoxLayout(group_box)
        group_layout.setSpacing(12)
        group_layout.setContentsMargins(16, 24, 16, 16)

        # Site Energy Calculation
        self.site_energy_button = QPushButton("Calculate Site Energies")
        self.site_energy_button.setMinimumHeight(40)
        self.site_energy_button.setProperty("class", "primary")
        group_layout.addWidget(self.site_energy_button)

        # CDC Analysis
        self.cdc_analysis_button = QPushButton("Run CDC Analysis")
        self.cdc_analysis_button.setMinimumHeight(40)
        self.cdc_analysis_button.setProperty("class", "secondary")
        group_layout.addWidget(self.cdc_analysis_button)

        # Info/Status area
        self.info_label = QLabel("💡 Load a PDB file to begin calculations")
        self.info_label.setStyleSheet("""
            QLabel {
                color: black !important;
                font-size: 13px;
                font-style: italic;
                margin-top: 10px;
                padding: 10px;
                background-color: #f3f3f3;
                border: 1px solid #e1e1e1;
                border-radius: 2px;
            }
        """)
        self.info_label.setWordWrap(True)
        group_layout.addWidget(self.info_label)

        layout.addWidget(group_box)
        
        # Progress information
        progress_group = QGroupBox("📊 Calculation Progress")
        progress_layout = QVBoxLayout(progress_group)
        
        self.progress_label = QLabel("No calculations running")
        self.progress_label.setStyleSheet("""
            QLabel {
                color: black !important;
                font-size: 12px;
                padding: 8px;
                background-color: #f3f3f3;
                border-radius: 2px;
            }
        """)
        progress_layout.addWidget(self.progress_label)
        
        layout.addWidget(progress_group)
        
        layout.addStretch()

    def connect_signals(self):
        """Connect widget signals"""
        self.site_energy_button.clicked.connect(self.calculate_site_energies.emit)
        self.cdc_analysis_button.clicked.connect(self.run_cdc_analysis.emit)

    def update_button_states(self):
        """Update button enabled/disabled states based on internal state"""
        self.site_energy_button.setEnabled(self.system_loaded)
        self.cdc_analysis_button.setEnabled(self.system_loaded and self.site_energies_calculated)

    def set_system_loaded(self, loaded: bool):
        """Set the system loaded state"""
        self.system_loaded = loaded
        self.update_button_states()
        
        if loaded:
            self.info_label.setText("✅ Structure loaded successfully\n💡 Click 'Calculate Site Energies' to begin")
        else:
            self.info_label.setText("💡 Load a PDB file to begin calculations")

    def set_site_energies_calculated(self, calculated: bool):
        """Set the site energies calculated state"""
        self.site_energies_calculated = calculated
        self.update_button_states()
        
        if calculated:
            self.info_label.setText("✅ Site energies calculated successfully\n🎯 Click 'Run CDC Analysis' for detailed contribution analysis")
        elif self.system_loaded:
            self.info_label.setText("✅ Structure loaded successfully\n💡 Click 'Calculate Site Energies' to begin")

    def set_calculation_running(self, running: bool, calc_type: str = ""):
        """Update UI during calculation"""
        if running:
            self.site_energy_button.setEnabled(False)
            self.cdc_analysis_button.setEnabled(False)
            if "site" in calc_type.lower():
                self.site_energy_button.setText("Calculating...")
                self.progress_label.setText("🔄 Calculating site energies...")
            elif "cdc" in calc_type.lower():
                self.cdc_analysis_button.setText("Running Analysis...")
                self.progress_label.setText("🔄 Running CDC analysis...")
        else:
            self.site_energy_button.setText("Calculate Site Energies")
            self.cdc_analysis_button.setText("Run CDC Analysis")
            self.progress_label.setText("No calculations running")
            self.update_button_states()
    
    def update_progress(self, percentage: int, message: str):
        """Update progress information"""
        self.progress_label.setText(f"🔄 {message} ({percentage}%)")
