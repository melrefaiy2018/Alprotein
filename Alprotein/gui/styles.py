"""
Backwards-compatible shim around :mod:`Alprotein.gui.theme`.

The constants below are kept so existing widgets that import them keep working.
New code should import the :class:`~Alprotein.gui.theme.Theme` objects directly.
"""

from __future__ import annotations

from .theme import light_theme

# Legacy token names — resolve to the light theme.
FONT_FAMILY = light_theme.font_family
FONT_SIZE_BODY = light_theme.font_size_body
FONT_SIZE_TITLE_PAGE = light_theme.font_size_page_title
FONT_SIZE_TITLE_CARD = light_theme.font_size_card_title
FONT_SIZE_LABEL = light_theme.font_size_label
FONT_SIZE_AXIS = "13px"

SPACING_UNIT = 4
PADDING_CARD = light_theme.padding_card
GAP_GRID = "16px"
SPACING_SECTION = "24px"

COLOR_PRIMARY = light_theme.text_primary
COLOR_SECONDARY = light_theme.text_secondary
COLOR_BG_APP = light_theme.bg_app
COLOR_BG_PANEL = light_theme.bg_panel
COLOR_BORDER = light_theme.border
COLOR_TEXT_PRIMARY = light_theme.text_primary
COLOR_TEXT_SECONDARY = light_theme.text_secondary
COLOR_SUCCESS = light_theme.success

# Composite outputs
GLOBAL_STYLESHEET = light_theme.qss()
CHART_STYLE = light_theme.chart_style()
