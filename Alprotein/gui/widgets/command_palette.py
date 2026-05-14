"""
A lightweight ⌘K-style command palette.

Register commands once at app startup, then trigger the palette with a global
shortcut. Matching is a simple subsequence fuzzy filter — good enough for a few
dozen actions and avoids any extra dependency.

    palette = CommandPalette(main_window)
    palette.add_command("file.open", "Open PDB…", main_window.on_open, shortcut="⌘O")
    palette.add_command("calc.spectrum", "Calculate absorption spectrum",
                        main_window.run_spectrum, group="Calculate")

    # Bound to Ctrl+K / ⌘K on the main window:
    QShortcut(QKeySequence("Ctrl+K"), main_window, palette.toggle)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QKeySequence
from PyQt5.QtWidgets import (
    QDialog,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QShortcut,
    QVBoxLayout,
    QWidget,
)

from ..theme import Theme, light_theme


@dataclass
class Command:
    """A palette entry."""

    key: str
    title: str
    callback: Callable[[], None]
    group: str = ""
    shortcut: str = ""
    keywords: List[str] = field(default_factory=list)

    def haystack(self) -> str:
        return " ".join([self.title, self.group, *self.keywords]).lower()


def _fuzzy_match(query: str, haystack: str) -> Optional[int]:
    """Return a score (lower = better) if ``query`` matches ``haystack`` as a subsequence."""
    if not query:
        return 0
    q = query.lower()
    h = haystack
    i = 0
    score = 0
    last = -1
    for ch in q:
        idx = h.find(ch, i)
        if idx < 0:
            return None
        if last >= 0:
            score += idx - last  # prefer adjacent characters
        last = idx
        i = idx + 1
    return score


class CommandPalette(QDialog):
    """Fuzzy-search dialog over registered actions."""

    def __init__(self, parent: QWidget, theme: Theme = light_theme) -> None:
        super().__init__(parent)
        self._commands: Dict[str, Command] = {}
        self._theme = theme

        self.setWindowFlags(Qt.Popup | Qt.FramelessWindowHint)
        self.setModal(True)
        self.resize(560, 360)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        self._search = QLineEdit()
        self._search.setPlaceholderText("Type a command…   (↑/↓ to navigate, Enter to run, Esc to close)")
        self._search.setClearButtonEnabled(True)
        self._search.textChanged.connect(self._refresh)
        layout.addWidget(self._search)

        self._list = QListWidget()
        self._list.itemActivated.connect(self._run_current)
        layout.addWidget(self._list, 1)

        self._apply_styles()

        # Esc closes; Enter on search bar runs.
        QShortcut(QKeySequence("Escape"), self, self.close)
        self._search.returnPressed.connect(self._run_current)

    # Public API -----------------------------------------------------------

    def add_command(
        self,
        key: str,
        title: str,
        callback: Callable[[], None],
        *,
        group: str = "",
        shortcut: str = "",
        keywords: Optional[List[str]] = None,
    ) -> None:
        self._commands[key] = Command(
            key=key,
            title=title,
            callback=callback,
            group=group,
            shortcut=shortcut,
            keywords=keywords or [],
        )

    def remove_command(self, key: str) -> None:
        self._commands.pop(key, None)

    def toggle(self) -> None:
        if self.isVisible():
            self.close()
        else:
            self.open_centered()

    def open_centered(self) -> None:
        parent = self.parentWidget()
        if parent is not None:
            geo = parent.geometry()
            x = geo.x() + (geo.width() - self.width()) // 2
            y = geo.y() + max(80, (geo.height() - self.height()) // 4)
            self.move(x, y)
        self._search.clear()
        self._refresh("")
        self.show()
        self._search.setFocus()

    # Internal -------------------------------------------------------------

    def _apply_styles(self) -> None:
        t = self._theme
        self.setStyleSheet(
            f"""
            QDialog {{
                background-color: {t.bg_panel};
                border: 1px solid {t.border_strong};
                border-radius: 10px;
            }}
            QLineEdit {{
                border: none;
                border-bottom: 1px solid {t.border};
                padding: 14px 16px;
                font-size: 14px;
                background-color: {t.bg_panel};
                color: {t.text_primary};
            }}
            QListWidget {{
                background-color: {t.bg_panel};
                border: none;
                padding: 6px;
                outline: none;
            }}
            QListWidget::item {{
                padding: 8px 10px;
                border-radius: 6px;
                color: {t.text_primary};
            }}
            QListWidget::item:selected {{
                background-color: {t.accent_soft};
                color: {t.text_primary};
            }}
            """
        )

    def _refresh(self, query: str = "") -> None:
        query = query or self._search.text()
        results: List[tuple[int, Command]] = []
        for cmd in self._commands.values():
            score = _fuzzy_match(query, cmd.haystack())
            if score is not None:
                results.append((score, cmd))
        results.sort(key=lambda pair: (pair[0], pair[1].title.lower()))

        self._list.clear()
        for _, cmd in results[:60]:
            text = cmd.title
            if cmd.group:
                text = f"{cmd.group} ›  {cmd.title}"
            if cmd.shortcut:
                text = f"{text}    {cmd.shortcut}"
            item = QListWidgetItem(text)
            item.setData(Qt.UserRole, cmd.key)
            self._list.addItem(item)
        if self._list.count() > 0:
            self._list.setCurrentRow(0)

    def _run_current(self, *_args) -> None:
        item = self._list.currentItem()
        if item is None and self._list.count() > 0:
            item = self._list.item(0)
        if item is None:
            return
        key = item.data(Qt.UserRole)
        cmd = self._commands.get(key)
        self.close()
        if cmd is not None:
            cmd.callback()

    # Key navigation -------------------------------------------------------

    def keyPressEvent(self, event):  # noqa: N802 — Qt signature
        if event.key() in (Qt.Key_Down, Qt.Key_Up):
            row = self._list.currentRow()
            if event.key() == Qt.Key_Down:
                row = min(row + 1, self._list.count() - 1)
            else:
                row = max(row - 1, 0)
            self._list.setCurrentRow(row)
            event.accept()
            return
        super().keyPressEvent(event)


__all__ = ["Command", "CommandPalette"]
