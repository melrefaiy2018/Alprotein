"""
Shared helpers for building the Alprotein workbench's card UI.

A "card" is a rounded, bordered panel with an optional title row at the top.
Qt's ``QGroupBox`` titles sit awkwardly *on top of* the border on macOS
(Big Sur and later) — they float in the widget's margin instead of cleanly
overlapping the border. That's the "PROJECT label floating outside the
card" effect.

Using ``QFrame`` + an internal ``QLabel`` instead gives us full control over
the layout, no platform-specific quirks, and matches the HTML mockup
1-for-1. All workbench panels migrate onto this helper.
"""

from __future__ import annotations

from typing import Optional, Tuple

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QFrame, QLabel, QVBoxLayout, QWidget


def make_card(
    title: str = "",
    *,
    parent: Optional[QWidget] = None,
    margins: Tuple[int, int, int, int] = (14, 12, 14, 14),
    spacing: int = 8,
) -> Tuple[QFrame, QVBoxLayout]:
    """Build a card frame with an optional title row.

    Returns ``(card, body_layout)`` so callers can simply
    ``body_layout.addWidget(...)`` their own content. The card itself
    already has ``class="card"`` set so the theme stylesheet styles it.
    """
    card = QFrame(parent)
    card.setProperty("class", "card")
    card.setFrameShape(QFrame.NoFrame)

    outer = QVBoxLayout(card)
    outer.setContentsMargins(*margins)
    outer.setSpacing(spacing)

    if title:
        label = QLabel(title)
        label.setProperty("class", "card-title")
        label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        outer.addWidget(label)

    return card, outer


__all__ = ["make_card"]
