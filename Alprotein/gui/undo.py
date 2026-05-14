"""
Undo / Redo support for the workbench.

The model is *input-level* undo: parameter edits, site-energy overrides,
axis-range changes — all reversible. Computed artefacts (Hamiltonian,
spectra) are not on the stack because rerunning them is the only safe way
to keep state consistent with inputs.

Usage::

    self.undo = UndoController(self)
    self.undo.push(ParameterEditCommand(self, key="dielectric_cdc",
                                        old=2.0, new=2.5))
"""

from __future__ import annotations

from typing import Any, Callable, Dict, Optional

from PyQt5.QtCore import QObject
from PyQt5.QtWidgets import QAction, QUndoCommand, QUndoStack, QWidget


# ---------------------------------------------------------------------------
# Concrete commands
# ---------------------------------------------------------------------------

class ParameterEditCommand(QUndoCommand):
    """Undoable change of a single calculation parameter."""

    def __init__(
        self,
        apply_fn: Callable[[str, Any], None],
        key: str,
        old: Any,
        new: Any,
        label: Optional[str] = None,
    ) -> None:
        super().__init__(label or f"Edit {key}")
        self._apply = apply_fn
        self._key = key
        self._old = old
        self._new = new

    def undo(self) -> None:  # noqa: D401 — Qt method
        self._apply(self._key, self._old)

    def redo(self) -> None:  # noqa: D401 — Qt method
        self._apply(self._key, self._new)

    def id(self) -> int:  # noqa: D401 — Qt method
        return hash(("parameter", self._key)) & 0x7FFFFFFF

    def mergeWith(self, other: QUndoCommand) -> bool:  # noqa: N802 — Qt method
        """Coalesce rapid edits to the same key (slider drag → one undo step)."""
        if isinstance(other, ParameterEditCommand) and other._key == self._key:
            self._new = other._new
            return True
        return False


class SiteEnergyOverrideCommand(QUndoCommand):
    """Undoable change to a single site-energy override."""

    SENTINEL = object()

    def __init__(
        self,
        apply_fn: Callable[[str, Optional[float]], None],
        pigment_id: str,
        old: Optional[float],
        new: Optional[float],
    ) -> None:
        super().__init__(f"Override {pigment_id}")
        self._apply = apply_fn
        self._pigment_id = pigment_id
        self._old = old
        self._new = new

    def undo(self) -> None:
        self._apply(self._pigment_id, self._old)

    def redo(self) -> None:
        self._apply(self._pigment_id, self._new)

    def id(self) -> int:
        return hash(("override", self._pigment_id)) & 0x7FFFFFFF

    def mergeWith(self, other: QUndoCommand) -> bool:  # noqa: N802 — Qt method
        if isinstance(other, SiteEnergyOverrideCommand) and other._pigment_id == self._pigment_id:
            self._new = other._new
            return True
        return False


# ---------------------------------------------------------------------------
# Controller
# ---------------------------------------------------------------------------

class UndoController(QObject):
    """Owns the :class:`QUndoStack` and exposes Edit-menu actions."""

    def __init__(self, parent: QWidget, undo_limit: int = 100) -> None:
        super().__init__(parent)
        self._stack = QUndoStack(self)
        self._stack.setUndoLimit(undo_limit)

    # Actions --------------------------------------------------------------

    def make_undo_action(self, parent: QWidget) -> QAction:
        action = self._stack.createUndoAction(parent, "Undo")
        action.setShortcut("Ctrl+Z")
        return action

    def make_redo_action(self, parent: QWidget) -> QAction:
        action = self._stack.createRedoAction(parent, "Redo")
        action.setShortcut("Ctrl+Shift+Z")
        return action

    # Stack ----------------------------------------------------------------

    def push(self, command: QUndoCommand) -> None:
        self._stack.push(command)

    def clear(self) -> None:
        self._stack.clear()

    @property
    def stack(self) -> QUndoStack:
        return self._stack


__all__ = [
    "ParameterEditCommand",
    "SiteEnergyOverrideCommand",
    "UndoController",
]
