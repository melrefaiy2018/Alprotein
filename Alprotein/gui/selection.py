"""
Central selection state shared across the workbench.

Every panel that needs to know "which pigment / exciton is the user looking
at?" subscribes to :class:`SelectionModel` instead of holding its own copy.
The model is mutated through ``set_pigment`` / ``set_exciton`` and emits
signals so dependent widgets can refresh.

Design notes
------------

* ``None`` means "nothing selected" — explicitly distinct from any pigment id.
* Setting the same value twice is a no-op (no signal emitted).
* The source of the change is carried in the signal so a widget can avoid
  reacting to its own emissions (e.g. the Hamiltonian heatmap shouldn't
  reselect itself when it just told the model about a click).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from PyQt5.QtCore import QObject, pyqtSignal


@dataclass(frozen=True)
class Selection:
    """Immutable snapshot of the current selection state."""

    pigment_id: Optional[str] = None
    exciton_k: Optional[int] = None  # 1-based index of an exciton eigenstate


class SelectionModel(QObject):
    """Observable selection state."""

    # source: short string identifying who initiated the change ("3d",
    # "hamiltonian", "spectrum", "inspector", "api"). Use it to short-circuit
    # feedback loops in the receiving widget.
    pigment_changed = pyqtSignal(object, str)   # pigment_id_or_None, source
    exciton_changed = pyqtSignal(object, str)   # k_or_None, source
    changed = pyqtSignal(object, str)           # Selection, source

    def __init__(self, parent: Optional[QObject] = None) -> None:
        super().__init__(parent)
        self._current = Selection()

    # ------------------------------------------------------------------
    # Read
    # ------------------------------------------------------------------

    @property
    def current(self) -> Selection:
        return self._current

    @property
    def pigment_id(self) -> Optional[str]:
        return self._current.pigment_id

    @property
    def exciton_k(self) -> Optional[int]:
        return self._current.exciton_k

    # ------------------------------------------------------------------
    # Mutate
    # ------------------------------------------------------------------

    def set_pigment(self, pigment_id: Optional[str], source: str = "api") -> None:
        if pigment_id == self._current.pigment_id:
            return
        self._current = Selection(pigment_id=pigment_id, exciton_k=self._current.exciton_k)
        self.pigment_changed.emit(pigment_id, source)
        self.changed.emit(self._current, source)

    def set_exciton(self, k: Optional[int], source: str = "api") -> None:
        if k == self._current.exciton_k:
            return
        self._current = Selection(pigment_id=self._current.pigment_id, exciton_k=k)
        self.exciton_changed.emit(k, source)
        self.changed.emit(self._current, source)

    def clear(self, source: str = "api") -> None:
        if self._current == Selection():
            return
        self._current = Selection()
        self.pigment_changed.emit(None, source)
        self.exciton_changed.emit(None, source)
        self.changed.emit(self._current, source)


__all__ = ["Selection", "SelectionModel"]
