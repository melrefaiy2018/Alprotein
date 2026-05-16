"""
Shared helpers for building Alprotein workbench sidebar sections.

A section is a lightly separated panel with an optional title row at the top.
Qt's ``QGroupBox`` titles sit awkwardly *on top of* the border on macOS
(Big Sur and later) — they float in the widget's margin instead of cleanly
overlapping the border. That's the "PROJECT label floating outside the
card" effect.

Using ``QFrame`` + an internal ``QLabel`` gives us full control over layout
without platform-specific QGroupBox title quirks.
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
    """Build a sidebar section with an optional title row.

    Returns ``(frame, body_layout)`` so callers can simply
    ``body_layout.addWidget(...)`` their own content.
    """
    card = QFrame(parent)
    card.setProperty("class", "section")
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
