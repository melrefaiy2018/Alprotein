"""
Right-side Inspector panel.

Shows everything we know about the currently-selected pigment: residue,
site energy (calculated, override, total), transition dipole, top exciton
contribution, strongest coupling partner. Driven by
:class:`~Alprotein.gui.selection.SelectionModel`.

The panel does not own any data — it pulls from a host that exposes
:meth:`InspectorHost.snapshot_for_pigment`. This keeps the widget testable
and avoids dragging the entire workbench into a unit test.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Protocol

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import (
    QDoubleSpinBox,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)


@dataclass(frozen=True)
class PigmentSnapshot:
    """Read-only view of one pigment's state, supplied by the host."""

    pigment_id: str
    resname: str = ""
    chain: str = ""
    resi: str = ""
    site_energy_cm: Optional[float] = None
    calculated_energy_cm: Optional[float] = None
    override_energy_cm: Optional[float] = None
    dipole_magnitude: Optional[float] = None
    domain: Optional[int] = None
    top_exciton: Optional[int] = None
    top_exciton_weight: Optional[float] = None
    top_partner_id: Optional[str] = None
    top_partner_coupling: Optional[float] = None


class InspectorHost(Protocol):
    """Anything that can serve pigment snapshots to the Inspector."""

    def snapshot_for_pigment(self, pigment_id: str) -> Optional[PigmentSnapshot]: ...


def _fmt(value: Optional[float], suffix: str = "", digits: int = 1) -> str:
    if value is None:
        return "—"
    return f"{value:.{digits}f}{suffix}"


class InspectorPanel(QWidget):
    """Sidebar that shows details for the currently-selected pigment."""

    focus_requested = pyqtSignal(str)               # pigment_id — focus camera
    override_changed = pyqtSignal(str, object)      # pigment_id, new_value | None

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self._host: Optional[InspectorHost] = None
        self._current_id: Optional[str] = None
        self._setup_ui()
        self.set_pigment(None)

    # ------------------------------------------------------------------
    # Host wiring
    # ------------------------------------------------------------------

    def set_host(self, host: InspectorHost) -> None:
        self._host = host

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _setup_ui(self) -> None:
        self.setProperty("class", "sidebar")
        outer = QVBoxLayout(self)
        outer.setContentsMargins(12, 12, 12, 12)
        outer.setSpacing(12)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)
        outer.addWidget(scroll, 1)

        body = QWidget()
        scroll.setWidget(body)
        v = QVBoxLayout(body)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(12)

        # --- Header card -------------------------------------------------
        header = QGroupBox("🧭 INSPECTOR")
        header.setProperty("class", "card")
        h_layout = QVBoxLayout(header)
        h_layout.setSpacing(6)
        self.title_label = QLabel("No selection")
        self.title_label.setStyleSheet("font-weight: 600; font-size: 14px;")
        h_layout.addWidget(self.title_label)

        self.subtitle_label = QLabel("Click a pigment in the 3D viewer or Hamiltonian.")
        self.subtitle_label.setProperty("class", "label")
        self.subtitle_label.setWordWrap(True)
        h_layout.addWidget(self.subtitle_label)

        self.focus_btn = QPushButton("Focus camera on pigment")
        self.focus_btn.setProperty("class", "primary")
        self.focus_btn.setEnabled(False)
        self.focus_btn.clicked.connect(self._emit_focus)
        h_layout.addWidget(self.focus_btn)

        v.addWidget(header)

        # --- Identity / structure ---------------------------------------
        v.addWidget(self._build_kv_card("Structure", [
            ("Chain / residue", "chain_resi"),
            ("Residue type", "resname"),
            ("Domain", "domain"),
        ]))

        # --- Energies ---------------------------------------------------
        v.addWidget(self._build_kv_card("Energies (cm⁻¹)", [
            ("Calculated E", "calculated_energy"),
            ("Override ΔE", "override_delta"),
            ("Effective E", "effective_energy"),
            ("|μ| (Debye)", "dipole"),
        ]))

        # --- Exciton ----------------------------------------------------
        v.addWidget(self._build_kv_card("Excitonic", [
            ("Top exciton k", "top_exciton"),
            ("Weight |c|²", "top_weight"),
            ("Strongest partner", "top_partner"),
            ("Partner J (cm⁻¹)", "top_partner_j"),
        ]))

        # --- Override editor --------------------------------------------
        override_card = QGroupBox("Site-energy override")
        override_card.setProperty("class", "card")
        ov = QVBoxLayout(override_card)
        ov.setSpacing(6)

        self.override_spin = QDoubleSpinBox()
        self.override_spin.setRange(8_000.0, 22_000.0)
        self.override_spin.setDecimals(1)
        self.override_spin.setSuffix(" cm⁻¹")
        self.override_spin.setEnabled(False)
        ov.addWidget(self.override_spin)

        row = QHBoxLayout()
        self.apply_override_btn = QPushButton("Apply")
        self.apply_override_btn.setProperty("class", "primary")
        self.apply_override_btn.setEnabled(False)
        self.apply_override_btn.clicked.connect(self._emit_apply_override)
        row.addWidget(self.apply_override_btn)

        self.clear_override_btn = QPushButton("Reset")
        self.clear_override_btn.setEnabled(False)
        self.clear_override_btn.clicked.connect(self._emit_clear_override)
        row.addWidget(self.clear_override_btn)
        ov.addLayout(row)

        v.addWidget(override_card)
        v.addStretch(1)

    def _build_kv_card(self, title: str, rows: list[tuple[str, str]]) -> QGroupBox:
        card = QGroupBox(title)
        card.setProperty("class", "card")
        layout = QVBoxLayout(card)
        layout.setSpacing(4)
        for label_text, key in rows:
            row = QHBoxLayout()
            label = QLabel(label_text)
            label.setProperty("class", "label")
            value = QLabel("—")
            value.setStyleSheet("font-weight: 600;")
            value.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
            value.setObjectName(f"inspector_value_{key}")
            row.addWidget(label, 1)
            row.addWidget(value, 1, Qt.AlignRight)
            layout.addLayout(row)
        return card

    # ------------------------------------------------------------------
    # Public slots
    # ------------------------------------------------------------------

    def set_pigment(self, pigment_id: Optional[str]) -> None:
        """Refresh the panel for the given pigment."""
        self._current_id = pigment_id
        snap = None
        if pigment_id and self._host is not None:
            snap = self._host.snapshot_for_pigment(pigment_id)

        if snap is None:
            self.title_label.setText("No selection")
            self.subtitle_label.setText("Click a pigment in the 3D viewer or Hamiltonian.")
            self._set_value("chain_resi", "—")
            self._set_value("resname", "—")
            self._set_value("domain", "—")
            self._set_value("calculated_energy", "—")
            self._set_value("override_delta", "—")
            self._set_value("effective_energy", "—")
            self._set_value("dipole", "—")
            self._set_value("top_exciton", "—")
            self._set_value("top_weight", "—")
            self._set_value("top_partner", "—")
            self._set_value("top_partner_j", "—")
            self.focus_btn.setEnabled(False)
            self.override_spin.setEnabled(False)
            self.apply_override_btn.setEnabled(False)
            self.clear_override_btn.setEnabled(False)
            return

        self.title_label.setText(snap.pigment_id)
        self.subtitle_label.setText(snap.resname or "Pigment")
        self._set_value("chain_resi", f"{snap.chain or '—'} / {snap.resi or '—'}")
        self._set_value("resname", snap.resname or "—")
        self._set_value("domain", "—" if snap.domain is None else f"D{snap.domain}")
        self._set_value("calculated_energy", _fmt(snap.calculated_energy_cm, "", 1))
        if snap.override_energy_cm is not None and snap.calculated_energy_cm is not None:
            delta = snap.override_energy_cm - snap.calculated_energy_cm
            self._set_value("override_delta", _fmt(delta, "", 1))
        else:
            self._set_value("override_delta", "—")
        self._set_value("effective_energy", _fmt(snap.site_energy_cm, "", 1))
        self._set_value("dipole", _fmt(snap.dipole_magnitude, " D", 2))
        self._set_value("top_exciton", "—" if snap.top_exciton is None else f"k = {snap.top_exciton}")
        self._set_value("top_weight", _fmt(snap.top_exciton_weight, "", 3))
        self._set_value("top_partner", snap.top_partner_id or "—")
        self._set_value("top_partner_j", _fmt(snap.top_partner_coupling, "", 1))

        self.focus_btn.setEnabled(True)
        self.override_spin.setEnabled(True)
        if snap.site_energy_cm is not None:
            self.override_spin.blockSignals(True)
            self.override_spin.setValue(float(snap.site_energy_cm))
            self.override_spin.blockSignals(False)
        self.apply_override_btn.setEnabled(snap.calculated_energy_cm is not None)
        self.clear_override_btn.setEnabled(snap.override_energy_cm is not None)

    def refresh(self) -> None:
        """Re-fetch the snapshot for the currently-shown pigment."""
        self.set_pigment(self._current_id)

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _set_value(self, key: str, text: str) -> None:
        label = self.findChild(QLabel, f"inspector_value_{key}")
        if label is not None:
            label.setText(text)

    def _emit_focus(self) -> None:
        if self._current_id:
            self.focus_requested.emit(self._current_id)

    def _emit_apply_override(self) -> None:
        if self._current_id:
            self.override_changed.emit(self._current_id, float(self.override_spin.value()))

    def _emit_clear_override(self) -> None:
        if self._current_id:
            self.override_changed.emit(self._current_id, None)


__all__ = ["InspectorPanel", "PigmentSnapshot", "InspectorHost"]
