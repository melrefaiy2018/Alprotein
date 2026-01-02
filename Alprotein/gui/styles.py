"""
Global Design System for Alprotein Application
"""

# Typography
FONT_FAMILY = "Segoe UI, Roboto, Helvetica Neue, Arial, sans-serif"
FONT_SIZE_BODY = "14px"
FONT_SIZE_TITLE_PAGE = "22px"
FONT_SIZE_TITLE_CARD = "16px"
FONT_SIZE_LABEL = "12px"
FONT_SIZE_AXIS = "13px"

# Spacing
SPACING_UNIT = 4
PADDING_CARD = "16px"
GAP_GRID = "16px"
SPACING_SECTION = "24px"

# Colors
COLOR_PRIMARY = "#000000"
COLOR_SECONDARY = "#3a3a3a"
COLOR_BG_APP = "#f7f7f7"
COLOR_BG_PANEL = "#ffffff"
COLOR_BORDER = "#e1e1e1"
COLOR_TEXT_PRIMARY = "#111111"
COLOR_TEXT_SECONDARY = "#6b6b6b"
COLOR_SUCCESS = "#1f7a44"

# Global QSS Stylesheet
GLOBAL_STYLESHEET = f"""
    /* Global Font */
    QWidget {{
        font-family: {FONT_FAMILY};
        font-size: {FONT_SIZE_BODY};
        color: {COLOR_TEXT_PRIMARY};
    }}

    /* Main Window & Backgrounds */
    QMainWindow, QWidget#centralWidget {{
        background-color: {COLOR_BG_APP};
    }}

    /* Dialogs */
    QDialog, QMessageBox {{
        background-color: #ffffff;
    }}
    QDialog QLabel, QMessageBox QLabel {{
        color: {COLOR_TEXT_PRIMARY};
    }}

    /* Sidebar */
    .sidebar {{
        background-color: {COLOR_BG_APP};
    }}

    /* ScrollArea */
    QScrollArea {{
        background-color: transparent;
        border: none;
    }}
    QScrollArea > QWidget > QWidget {{
        background-color: transparent;
    }}
    /* Cards (QFrame/QGroupBox) */
    QFrame.card, QGroupBox.card {{
        background-color: {COLOR_BG_PANEL};
        border: 1px solid {COLOR_BORDER};
        border-radius: 14px;
        padding: {PADDING_CARD};
        margin-top: 16px; /* Space for title */
    }}
    
    QGroupBox.card::title {{
        subcontrol-origin: margin;
        subcontrol-position: top left;
        left: 16px;
        top: 0px;
        padding: 0 4px;
        background-color: {COLOR_BG_PANEL};
        color: {COLOR_TEXT_PRIMARY};
        font-size: {FONT_SIZE_TITLE_CARD};
        font-weight: 600;
    }}
    
    /* Card Titles */
    QLabel.card-title {{
        font-size: {FONT_SIZE_TITLE_CARD};
        font-weight: 600; /* Semibold */
        color: {COLOR_TEXT_PRIMARY};
        padding-bottom: 8px;
    }}

    /* Page Titles */
    QLabel.page-title {{
        font-size: {FONT_SIZE_TITLE_PAGE};
        font-weight: 600;
        color: {COLOR_TEXT_PRIMARY};
        margin-bottom: {SPACING_SECTION};
    }}

    /* Labels */
    QLabel.label {{
        color: {COLOR_TEXT_SECONDARY};
        font-size: {FONT_SIZE_LABEL};
    }}
    
    QLabel.label-success {{
        color: {COLOR_SUCCESS};
        font-weight: 600;
        font-size: {FONT_SIZE_LABEL};
    }}

    /* Buttons */
    QPushButton {{
        height: 36px;
        border-radius: 6px;
        font-weight: 500;
        padding: 0 16px;
    }}
    
    QPushButton.primary {{
        background-color: {COLOR_PRIMARY};
        color: #ffffff;
        border: none;
    }}
    QPushButton.primary:hover {{
        background-color: #000000;
    }}

    QPushButton.secondary {{
        background-color: transparent;
        color: {COLOR_PRIMARY};
        border: 1px solid {COLOR_BORDER};
    }}
    QPushButton.secondary:hover {{
        background-color: #f5f5f5;
        border-color: {COLOR_PRIMARY};
    }}

    /* Small Action Buttons (Details, Export) */
    QPushButton.action-small {{
        height: 32px;
        font-size: 12px;
        padding: 0 12px;
        background-color: #ffffff;
        color: {COLOR_PRIMARY};
        border: 1px solid {COLOR_BORDER};
        border-radius: 4px;
    }}
    QPushButton.action-small:hover {{
        background-color: #f5f5f5;
    }}

    /* Inputs */
    QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {{
        height: 36px;
        border: 1px solid {COLOR_BORDER};
        border-radius: 6px;
        padding: 0 8px;
        background-color: #ffffff;
    }}
    
    /* Tabs */
    QTabWidget::pane {{
        border: none;
        background: transparent;
    }}
    QTabBar::tab {{
        background: transparent;
        color: {COLOR_TEXT_SECONDARY};
        padding: 8px 16px;
        font-size: 14px;
        font-weight: 500;
        border-bottom: 2px solid transparent;
        margin-right: 4px;
    }}
    QTabBar::tab:selected {{
        color: {COLOR_PRIMARY};
        border-bottom: 2px solid {COLOR_PRIMARY};
    }}
    QTabBar::tab:hover {{
        color: {COLOR_PRIMARY};
    }}
    /* Tables */
    QTableWidget {{
        background-color: #ffffff;
        alternate-background-color: #f9fafb;
        gridline-color: #e5e7eb;
        border: 1px solid {COLOR_BORDER};
        border-radius: 6px;
        color: {COLOR_TEXT_PRIMARY};
    }}
    QHeaderView {{
        background-color: #f3f4f6;
        border: none;
    }}
    QHeaderView::section {{
        background-color: #f3f4f6;
        color: {COLOR_TEXT_PRIMARY};
        padding: 8px;
        border: none;
        border-bottom: 1px solid {COLOR_BORDER};
        border-right: 1px solid {COLOR_BORDER};
        font-weight: 600;
        font-size: 12px;
    }}
    QTableCornerButton::section {{
        background-color: #f3f4f6;
        border: none;
        border-bottom: 1px solid {COLOR_BORDER};
        border-right: 1px solid {COLOR_BORDER};
    }}
    QTableWidget::item {{
        padding: 6px;
        color: {COLOR_TEXT_PRIMARY};
    }}
    QTableWidget::item:selected {{
        background-color: {COLOR_PRIMARY};
        color: #ffffff;
    }}"""

# Chart Styling Configuration
CHART_STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Segoe UI", "Roboto", "Helvetica Neue", "Arial"],
    "font.size": 12,
    "axes.labelsize": 13,
    "axes.labelweight": "semibold",
    "axes.titlesize": 14,
    "axes.titleweight": "semibold",
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linewidth": 0.5,
    "lines.linewidth": 2.0,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "figure.autolayout": True,
}
