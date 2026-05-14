"""
Non-blocking toast notifications anchored to the bottom-right of a parent window.

Replaces ``QMessageBox.information``/``warning`` for status feedback that should
not interrupt the user. Use :class:`ToastManager` from the main window:

    self.toasts = ToastManager(self)
    self.toasts.info("Hamiltonian built", "13×13 · max |J| = 92.4 cm⁻¹")
    self.toasts.error("TrEsp parameters missing", details=traceback_text)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

from PyQt5.QtCore import (
    QEvent,
    QObject,
    QPoint,
    QPropertyAnimation,
    QEasingCurve,
    Qt,
    QTimer,
)
from PyQt5.QtGui import QGuiApplication
from PyQt5.QtWidgets import (
    QApplication,
    QFrame,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from ..theme import Theme, light_theme

Severity = str  # "info" | "success" | "warning" | "error"


@dataclass(frozen=True)
class _SeverityStyle:
    icon: str
    fg: str
    border: str


def _style_for(theme: Theme, severity: Severity) -> _SeverityStyle:
    mapping = {
        "info": _SeverityStyle("ℹ", theme.text_primary, theme.border_strong),
        "success": _SeverityStyle("✓", theme.success, theme.success),
        "warning": _SeverityStyle("⚠", theme.warning, theme.warning),
        "error": _SeverityStyle("✕", theme.danger, theme.danger),
    }
    return mapping.get(severity, mapping["info"])


class Toast(QFrame):
    """A single toast card. Self-dismisses after ``duration_ms``."""

    GAP = 8
    DEFAULT_DURATION_MS = 6000

    def __init__(
        self,
        parent: QWidget,
        title: str,
        body: str = "",
        severity: Severity = "info",
        duration_ms: Optional[int] = None,
        details: Optional[str] = None,
        theme: Theme = light_theme,
    ) -> None:
        super().__init__(parent)
        self._theme = theme
        self._details = details

        self.setWindowFlags(Qt.SubWindow | Qt.FramelessWindowHint)
        self.setAttribute(Qt.WA_TransparentForMouseEvents, False)
        self.setAttribute(Qt.WA_StyledBackground, True)
        self.setFixedWidth(340)

        style = _style_for(theme, severity)
        self.setStyleSheet(
            f"""
            QFrame {{
                background-color: {theme.bg_panel};
                border: 1px solid {style.border};
                border-radius: 8px;
            }}
            QLabel.title {{
                color: {style.fg};
                font-weight: 600;
                font-size: 13px;
            }}
            QLabel.body {{
                color: {theme.text_secondary};
                font-size: 12px;
            }}
            QPushButton.close {{
                background: transparent;
                color: {theme.text_secondary};
                border: none;
                font-size: 14px;
                padding: 0;
            }}
            QPushButton.close:hover {{
                color: {theme.text_primary};
            }}
            QPushButton.details {{
                background: transparent;
                color: {theme.accent};
                border: none;
                padding: 0;
                font-size: 11px;
                text-align: left;
            }}
            """
        )

        outer = QVBoxLayout(self)
        outer.setContentsMargins(12, 10, 12, 10)
        outer.setSpacing(4)

        header = QHBoxLayout()
        header.setSpacing(8)
        title_lbl = QLabel(f"{style.icon}  {title}")
        title_lbl.setProperty("class", "title")
        title_lbl.setStyleSheet(f"color: {style.fg}; font-weight: 600;")
        header.addWidget(title_lbl, 1)
        close_btn = QPushButton("×")
        close_btn.setProperty("class", "close")
        close_btn.setFixedSize(18, 18)
        close_btn.setCursor(Qt.PointingHandCursor)
        close_btn.clicked.connect(self.dismiss)
        header.addWidget(close_btn)
        outer.addLayout(header)

        if body:
            body_lbl = QLabel(body)
            body_lbl.setProperty("class", "body")
            body_lbl.setStyleSheet(f"color: {theme.text_secondary};")
            body_lbl.setWordWrap(True)
            outer.addWidget(body_lbl)

        if details:
            details_btn = QPushButton("Copy details ›")
            details_btn.setProperty("class", "details")
            details_btn.setStyleSheet(f"color: {theme.accent};")
            details_btn.setCursor(Qt.PointingHandCursor)
            details_btn.clicked.connect(self._copy_details)
            outer.addWidget(details_btn, 0, Qt.AlignLeft)

        self._duration_ms = duration_ms or self.DEFAULT_DURATION_MS
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.timeout.connect(self.dismiss)

    def show_at(self, target_pos: QPoint) -> None:
        """Position the toast at ``target_pos`` (top-left in parent coords) and start the dismiss timer."""
        self.move(target_pos)
        self.show()
        self.adjustSize()
        if self._duration_ms > 0:
            self._timer.start(self._duration_ms)

    def dismiss(self) -> None:
        """Fade and hide. The owning manager will reflow remaining toasts."""
        self._timer.stop()
        anim = QPropertyAnimation(self, b"windowOpacity", self)
        anim.setDuration(180)
        anim.setStartValue(1.0)
        anim.setEndValue(0.0)
        anim.setEasingCurve(QEasingCurve.OutCubic)
        anim.finished.connect(self.deleteLater)
        anim.start()
        owner = getattr(self, "_owner", None)
        if owner is not None:
            owner._on_toast_dismissed(self)

    def _copy_details(self) -> None:
        if not self._details:
            return
        cb = QGuiApplication.clipboard() or QApplication.clipboard()
        if cb is not None:
            cb.setText(self._details)


class ToastManager:
    """Owns and stacks active toasts inside a parent widget."""

    MARGIN = 18

    def __init__(self, parent: QWidget, theme: Theme = light_theme) -> None:
        self._parent = parent
        self._theme = theme
        self._active: List[Toast] = []
        parent.installEventFilter(_ResizeRelayout(parent, self))

    def set_theme(self, theme: Theme) -> None:
        self._theme = theme

    # Public API -----------------------------------------------------------

    def info(self, title: str, body: str = "", **kwargs) -> Toast:
        return self._show(title, body, "info", **kwargs)

    def success(self, title: str, body: str = "", **kwargs) -> Toast:
        return self._show(title, body, "success", **kwargs)

    def warning(self, title: str, body: str = "", **kwargs) -> Toast:
        return self._show(title, body, "warning", **kwargs)

    def error(
        self,
        title: str,
        body: str = "",
        details: Optional[str] = None,
        **kwargs,
    ) -> Toast:
        # Errors stay longer and offer a copy-details affordance.
        kwargs.setdefault("duration_ms", 10_000)
        return self._show(title, body, "error", details=details, **kwargs)

    # Internal -------------------------------------------------------------

    def _show(self, title: str, body: str, severity: Severity, **kwargs) -> Toast:
        toast = Toast(self._parent, title, body, severity, theme=self._theme, **kwargs)
        # Replace owner so dismiss callback finds us, while keeping parent for coords.
        toast._owner = self  # type: ignore[attr-defined]
        self._active.append(toast)
        self._layout()
        toast.show_at(self._compute_pos(toast))
        return toast

    def _on_toast_dismissed(self, toast: Toast) -> None:
        if toast in self._active:
            self._active.remove(toast)
        self._layout()

    def _layout(self) -> None:
        """Reflow stacked toasts so they appear from bottom up."""
        for toast in list(self._active):
            if not toast.isVisible():
                continue
            toast.move(self._compute_pos(toast))

    def _compute_pos(self, toast: Toast) -> QPoint:
        idx = self._active.index(toast)
        parent_rect = self._parent.rect()
        x = parent_rect.right() - toast.width() - self.MARGIN
        y = parent_rect.bottom() - self.MARGIN
        for above in self._active[idx:]:
            y -= above.sizeHint().height() + Toast.GAP
        return QPoint(x, y)


class _ResizeRelayout(QObject):
    """Event filter that relayouts toasts when the parent resizes or moves."""

    def __init__(self, parent: QWidget, manager: ToastManager) -> None:
        super().__init__(parent)
        self._parent = parent
        self._manager = manager

    def eventFilter(self, obj, event):  # noqa: N802 — Qt signature
        if obj is self._parent and event.type() in (QEvent.Resize, QEvent.Move):
            self._manager._layout()
        return False


__all__ = ["Toast", "ToastManager"]
