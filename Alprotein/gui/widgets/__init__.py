"""GUI widgets for the Alprotein Scientific Workbench."""

from .analysis_dashboard import AnalysisDashboard
from .data_table_widget import DataTableWidget
from .hamiltonian_widget import HamiltonianWidget
from .protein_viewer import ProteinViewer
from .spectrum_widget import SpectrumPlotWidget
from .summary_cards_widget import SummaryCardsWidget
from .tools_panel import ToolsPanel

__all__ = [
    "AnalysisDashboard",
    "DataTableWidget",
    "HamiltonianWidget",
    "ProteinViewer",
    "SpectrumPlotWidget",
    "SummaryCardsWidget",
    "ToolsPanel",
]
