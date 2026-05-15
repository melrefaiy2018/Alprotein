"""
Theme tokens for the Alprotein Scientific Workbench.

Two themes are provided: :data:`light_theme` (default) and :data:`dark_theme`.
A :class:`Theme` is a frozen dataclass of tokens — colors, typography, spacing —
plus a :meth:`Theme.qss` method that renders a Qt stylesheet, and a
:meth:`Theme.chart_style` method that returns matplotlib rcParams.

The legacy module ``Alprotein.gui.styles`` re-exports :data:`GLOBAL_STYLESHEET`
and :data:`CHART_STYLE` for backwards compatibility — both resolve to the light
theme.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field, replace
from typing import Dict


def _platform_font_stack() -> str:
    """Return a font-family chain whose first hit is installed on this OS.

    Qt logs a warning and burns ~150 ms scanning aliases for every family
    name it can't resolve, so we keep the stack minimal and OS-specific.
    """
    if sys.platform == "darwin":
        # ``.AppleSystemUIFont`` resolves to SF UI on Big Sur+.
        return "'.AppleSystemUIFont', 'Helvetica Neue', Helvetica, Arial"
    if sys.platform.startswith("win"):
        return "'Segoe UI', Arial"
    return "'Noto Sans', 'DejaVu Sans', Arial"


@dataclass(frozen=True)
class Theme:
    """Visual tokens for one theme variant."""

    name: str

    # Background & surface
    bg_app: str
    bg_panel: str
    bg_subtle: str

    # Foreground
    text_primary: str
    text_secondary: str
    text_inverse: str

    # Borders
    border: str
    border_strong: str

    # Brand / accent
    accent: str
    accent_hover: str
    accent_soft: str

    # Semantic
    success: str
    success_soft: str
    warning: str
    warning_soft: str
    danger: str
    danger_soft: str

    # Typography — resolved at import time for the running OS so Qt
    # doesn't burn ~150 ms scanning aliases that aren't installed (e.g.
    # "Segoe UI" on macOS). See _platform_font_stack().
    font_family: str = field(default_factory=_platform_font_stack)
    font_size_body: str = "13px"
    font_size_label: str = "12px"
    font_size_card_title: str = "16px"
    font_size_page_title: str = "22px"

    # Spacing & radius
    radius_card: str = "10px"
    radius_input: str = "6px"
    radius_pill: str = "999px"
    padding_card: str = "14px"

    def qss(self) -> str:
        """Render the full QSS stylesheet for this theme."""
        t = self
        return f"""
        /* Global */
        QWidget {{
            font-family: {t.font_family};
            font-size: {t.font_size_body};
            color: {t.text_primary};
        }}

        QMainWindow, QDialog {{
            background-color: {t.bg_app};
        }}

        /* Sidebar surfaces */
        .sidebar, QWidget#sidebar {{
            background-color: {t.bg_subtle};
        }}

        /* Scrollarea */
        QScrollArea {{
            background-color: transparent;
            border: none;
        }}
        QScrollArea > QWidget > QWidget {{
            background-color: transparent;
        }}

        /* Cards — preferred style: QFrame.card with an internal
           QLabel.card-title. QGroupBox.card stays supported for the few
           legacy widgets that still use it, but its native title placement
           on macOS misbehaves (floats outside the border), so new code
           should use widgets.cards.make_card(). */
        QFrame.card, QGroupBox.card {{
            background-color: {t.bg_panel};
            border: 1px solid {t.border};
            border-radius: {t.radius_card};
        }}
        QGroupBox.card {{
            padding: {t.padding_card};
            margin-top: 16px;
        }}
        QGroupBox.card::title {{
            subcontrol-origin: margin;
            subcontrol-position: top left;
            left: 14px;
            top: -4px;
            padding: 0 6px;
            background-color: {t.bg_panel};
            color: {t.text_primary};
            font-size: {t.font_size_card_title};
            font-weight: 600;
        }}

        /* Card title — the label inserted by make_card() */
        QLabel.card-title {{
            font-size: 11px;
            font-weight: 700;
            letter-spacing: 1.5px;
            color: {t.text_secondary};
            text-transform: uppercase;
            padding-bottom: 4px;
        }}
        QLabel.page-title {{
            font-size: {t.font_size_page_title};
            font-weight: 600;
            color: {t.text_primary};
        }}
        QLabel.label {{
            color: {t.text_secondary};
            font-size: {t.font_size_label};
        }}
        QLabel.label-success {{
            color: {t.success};
            font-weight: 600;
            font-size: {t.font_size_label};
        }}

        /* Buttons */
        QPushButton {{
            height: 32px;
            border-radius: {t.radius_input};
            font-weight: 500;
            padding: 0 14px;
            background-color: {t.bg_panel};
            color: {t.text_primary};
            border: 1px solid {t.border_strong};
        }}
        QPushButton:hover {{
            background-color: {t.bg_subtle};
        }}
        QPushButton:disabled {{
            color: {t.text_secondary};
            background-color: {t.bg_subtle};
        }}
        QPushButton.primary {{
            background-color: {t.accent};
            color: {t.text_inverse};
            border: 1px solid {t.accent};
        }}
        QPushButton.primary:hover {{
            background-color: {t.accent_hover};
            border-color: {t.accent_hover};
        }}
        QPushButton.secondary {{
            background-color: transparent;
            color: {t.text_primary};
            border: 1px solid {t.border};
        }}
        QPushButton.secondary:hover {{
            background-color: {t.bg_subtle};
            border-color: {t.text_primary};
        }}
        QPushButton.danger {{
            color: {t.danger};
            border-color: {t.danger_soft};
        }}
        QPushButton.danger:hover {{
            background-color: {t.danger_soft};
        }}
        QPushButton.action-small {{
            height: 28px;
            font-size: 12px;
            padding: 0 10px;
        }}

        /* Inputs */
        QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {{
            height: 30px;
            border: 1px solid {t.border_strong};
            border-radius: {t.radius_input};
            padding: 0 8px;
            background-color: {t.bg_panel};
            color: {t.text_primary};
        }}
        QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {{
            border-color: {t.accent};
        }}

        /* Tabs */
        QTabWidget::pane {{
            border: none;
            background: transparent;
        }}
        QTabBar::tab {{
            background: transparent;
            color: {t.text_secondary};
            padding: 8px 14px;
            font-size: 13px;
            font-weight: 500;
            border-bottom: 2px solid transparent;
            margin-right: 2px;
        }}
        QTabBar::tab:selected {{
            color: {t.text_primary};
            border-bottom: 2px solid {t.accent};
        }}
        QTabBar::tab:hover:!selected {{
            color: {t.text_primary};
        }}

        /* Tables */
        QTableWidget {{
            background-color: {t.bg_panel};
            alternate-background-color: {t.bg_subtle};
            gridline-color: {t.border};
            border: 1px solid {t.border};
            border-radius: {t.radius_input};
            color: {t.text_primary};
        }}
        QHeaderView::section {{
            background-color: {t.bg_subtle};
            color: {t.text_primary};
            padding: 6px 8px;
            border: none;
            border-bottom: 1px solid {t.border};
            border-right: 1px solid {t.border};
            font-weight: 600;
            font-size: 12px;
        }}
        QTableWidget::item:selected {{
            background-color: {t.accent};
            color: {t.text_inverse};
        }}

        /* Progress */
        QProgressBar {{
            border: 1px solid {t.border};
            border-radius: 3px;
            text-align: center;
            background-color: {t.bg_subtle};
            color: {t.text_primary};
            height: 8px;
        }}
        QProgressBar::chunk {{
            background-color: {t.accent};
            border-radius: 3px;
        }}

        /* Tooltips */
        QToolTip {{
            background-color: {t.bg_panel};
            color: {t.text_primary};
            border: 1px solid {t.border_strong};
            padding: 6px 8px;
            border-radius: 4px;
        }}
        """

    def chart_style(self) -> Dict[str, object]:
        """Matplotlib rcParams that match this theme."""
        return {
            "font.family": "sans-serif",
            "font.sans-serif": ["Segoe UI", "Roboto", "Helvetica Neue", "Arial"],
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.labelweight": "semibold",
            "axes.titlesize": 13,
            "axes.titleweight": "semibold",
            "axes.edgecolor": self.border_strong,
            "axes.labelcolor": self.text_primary,
            "xtick.color": self.text_secondary,
            "ytick.color": self.text_secondary,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "axes.grid": True,
            "grid.alpha": 0.3,
            "grid.linewidth": 0.5,
            "grid.color": self.border,
            "lines.linewidth": 2.0,
            "figure.facecolor": self.bg_panel,
            "axes.facecolor": self.bg_panel,
            "figure.autolayout": True,
            "text.color": self.text_primary,
        }


# ---------------------------------------------------------------------------
# Variants
# ---------------------------------------------------------------------------

light_theme = Theme(
    name="light",
    bg_app="#f4f5f7",
    bg_panel="#ffffff",
    bg_subtle="#f7f8fa",
    text_primary="#14181f",
    text_secondary="#6b7280",
    text_inverse="#ffffff",
    border="#e3e6ea",
    border_strong="#c9ced6",
    accent="#2563eb",
    accent_hover="#1d4ed8",
    accent_soft="#dbeafe",
    success="#1f7a44",
    success_soft="#dcfce7",
    warning="#8a6d1f",
    warning_soft="#fef3c7",
    danger="#8b2b2b",
    danger_soft="#fee2e2",
)


dark_theme = replace(
    light_theme,
    name="dark",
    bg_app="#0f1115",
    bg_panel="#171a21",
    bg_subtle="#1f232c",
    text_primary="#e6e8eb",
    text_secondary="#9aa3af",
    text_inverse="#0f1115",
    border="#2a2f39",
    border_strong="#3a4150",
    accent="#3b82f6",
    accent_hover="#60a5fa",
    accent_soft="#1e3a8a",
    success="#34d399",
    success_soft="#064e3b",
    warning="#fbbf24",
    warning_soft="#78350f",
    danger="#f87171",
    danger_soft="#7f1d1d",
)


THEMES: Dict[str, Theme] = {
    light_theme.name: light_theme,
    dark_theme.name: dark_theme,
}


def get_theme(name: str) -> Theme:
    """Return a theme by name; falls back to the light theme."""
    return THEMES.get(name, light_theme)
