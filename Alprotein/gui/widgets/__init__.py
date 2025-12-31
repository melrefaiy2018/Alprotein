"""
GUI Widgets for Alprotein
"""

from .file_loader import FileLoaderWidget
from .calculation_panel import CalculationPanel
from .results_panel import ResultsPanel
from .protein_viewer import ProteinViewer
from .control_sidebar import ControlSidebar
from .actions_panel import ActionsPanel
from .control_panel_3d import ControlPanel3D

__all__ = [
    'FileLoaderWidget',
    'CalculationPanel', 
    'ResultsPanel',
    'ProteinViewer',
    'ControlSidebar',
    'ActionsPanel',
    'ControlPanel3D'
]
