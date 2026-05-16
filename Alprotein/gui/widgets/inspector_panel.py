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
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFrame,
    QHBoxLayout,
    QLabel,
    QPlainTextEdit,
    QPushButton,
    QScrollArea,
    QSlider,
    QVBoxLayout,
    QWidget,
)

from Alprotein.gui.widgets.cards import make_card


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
    notebook_edited = pyqtSignal(str)               # full notebook text
    color_mode_changed = pyqtSignal(str)            # "pigment_type" | "site_energy" | "exciton" | "domain"
    domain_cutoff_changed = pyqtSignal(int)         # cm⁻¹, integer-valued for slider
    domain_focus_requested = pyqtSignal(int)        # 1-based domain id

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
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        outer.addWidget(scroll, 1)

        body = QWidget()
        scroll.setWidget(body)
        v = QVBoxLayout(body)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(0)

        # --- Header card -------------------------------------------------
        header, h_layout = make_card("Inspector", margins=(18, 16, 18, 18), spacing=8)
        self.title_label = QLabel("No selection")
        self.title_label.setStyleSheet("font-weight: 700; font-size: 14px; color: #111827;")
        h_layout.addWidget(self.title_label)

        self.subtitle_label = QLabel("Click a pigment in the 3D viewer or Hamiltonian.")
        self.subtitle_label.setProperty("class", "label")
        self.subtitle_label.setWordWrap(True)
        h_layout.addWidget(self.subtitle_label)

        self.focus_btn = QPushButton("Focus camera")
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
        override_card, ov = make_card("Site-energy override", margins=(18, 16, 18, 18), spacing=8)

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

        # --- Live View Tools card ----------------------------------------
        view_tools_card, vt = make_card("Live view", margins=(18, 16, 18, 18), spacing=8)

        row = QHBoxLayout()
        row.addWidget(QLabel("Color by"))
        self._color_mode = QComboBox()
        self._color_mode.addItem("Pigment type", "pigment_type")
        self._color_mode.addItem("Site energy", "site_energy")
        self._color_mode.addItem("|c|² of exciton k", "exciton")
        self._color_mode.addItem("Domain", "domain")
        self._color_mode.currentIndexChanged.connect(
            lambda _i: self.color_mode_changed.emit(self._color_mode.currentData())
        )
        row.addWidget(self._color_mode, 1)
        vt.addLayout(row)

        self._show_dipoles = QCheckBox("Show μ vectors")
        self._show_dipoles.setEnabled(False)
        self._show_dipoles.setToolTip(
            "Render transition dipole vectors on each pigment.\n"
            "Coming with the dichroism work (P1.1)."
        )
        vt.addWidget(self._show_dipoles)

        self._show_labels = QCheckBox("Show pigment labels")
        self._show_labels.setChecked(True)
        self._show_labels.setEnabled(False)
        self._show_labels.setToolTip("Always on for now; toggle wiring is a TODO.")
        vt.addWidget(self._show_labels)

        v.addWidget(view_tools_card)

        # --- Domains card ------------------------------------------------
        self._domains_card, self._domains_layout = make_card("Domains", margins=(18, 16, 18, 18), spacing=8)
        # Track dynamically-rendered domain rows so set_domains() can clear
        # them without relying on layout indices.
        self._domain_row_widgets: list[QWidget] = []
        self._domains_empty = QLabel("Build the Hamiltonian to discover excitonic domains.")
        self._domains_empty.setProperty("class", "label")
        self._domains_empty.setWordWrap(True)
        self._domains_layout.addWidget(self._domains_empty)

        row = QHBoxLayout()
        row.addWidget(QLabel("Cutoff |J|"))
        self._domain_cutoff = QSlider(Qt.Horizontal)
        self._domain_cutoff.setRange(0, 100)
        self._domain_cutoff.setValue(20)
        self._domain_cutoff.valueChanged.connect(self.domain_cutoff_changed.emit)
        self._domain_cutoff_readout = QLabel("20 cm⁻¹")
        self._domain_cutoff.valueChanged.connect(
            lambda v: self._domain_cutoff_readout.setText(f"{v} cm⁻¹")
        )
        row.addWidget(self._domain_cutoff, 1)
        row.addWidget(self._domain_cutoff_readout)
        cutoff_widget = QWidget()
        cutoff_widget.setLayout(row)
        self._domains_cutoff_row = cutoff_widget
        self._domains_layout.addWidget(cutoff_widget)

        v.addWidget(self._domains_card)

        # --- Notebook card -----------------------------------------------
        notebook_card, nb = make_card("Notebook", margins=(18, 16, 18, 18), spacing=6)
        hint = QLabel("Per-project notes (saved inside .alproj).")
        hint.setProperty("class", "label")
        hint.setWordWrap(True)
        nb.addWidget(hint)

        self._notebook_edit = QPlainTextEdit()
        self._notebook_edit.setPlaceholderText("Type notes here…")
        self._notebook_edit.setMinimumHeight(120)
        self._notebook_edit.textChanged.connect(
            lambda: self.notebook_edited.emit(self._notebook_edit.toPlainText())
        )
        nb.addWidget(self._notebook_edit)

        v.addWidget(notebook_card)

        v.addStretch(1)

    # ------------------------------------------------------------------
    # Public setters used by the workbench
    # ------------------------------------------------------------------

    def set_notebook_text(self, text: str) -> None:
        """Programmatically set the notebook (e.g. after loading a project)."""
        if text == self._notebook_edit.toPlainText():
            return
        self._notebook_edit.blockSignals(True)
        self._notebook_edit.setPlainText(text or "")
        self._notebook_edit.blockSignals(False)

    def get_notebook_text(self) -> str:
        return self._notebook_edit.toPlainText()

    def set_domains(self, domains: dict, palette: Optional[dict] = None) -> None:
        """Render coloured swatches for the discovered domains.

        ``domains`` maps ``{domain_id: [pigment_index, ...]}``. If
        ``palette`` is supplied (``{domain_id: "#rrggbb"}``) we use it,
        otherwise we cycle a small built-in colour list.
        """
        # Drop the previous dynamic rows (we keep card title / empty hint /
        # cutoff slider intact via explicit tracking).
        for w in self._domain_row_widgets:
            self._domains_layout.removeWidget(w)
            w.deleteLater()
        self._domain_row_widgets.clear()

        if not domains:
            self._domains_empty.setVisible(True)
            return
        self._domains_empty.setVisible(False)

        default_palette = [
            "#2563eb", "#22c55e", "#f59e0b", "#ef4444",
            "#a855f7", "#06b6d4", "#84cc16", "#d946ef",
        ]
        # Insert each domain row just before the cutoff slider.
        cutoff_idx = self._domains_layout.indexOf(self._domains_cutoff_row)
        insert_at = cutoff_idx if cutoff_idx >= 0 else self._domains_layout.count()
        for i, (did, members) in enumerate(sorted(domains.items())):
            colour = (palette or {}).get(did, default_palette[i % len(default_palette)])
            row = QHBoxLayout()
            swatch = QLabel()
            swatch.setFixedSize(12, 12)
            swatch.setStyleSheet(
                f"background-color: {colour}; border-radius: 2px;"
            )
            name = QLabel(f"D{int(did) + 1}")
            name.setStyleSheet("font-weight: 600;")
            count = QLabel(f"{len(members)} pigments")
            count.setProperty("class", "label")
            row.addWidget(swatch)
            row.addWidget(name)
            row.addStretch(1)
            row.addWidget(count)
            container = QWidget()
            container.setLayout(row)
            container.mousePressEvent = (  # type: ignore[assignment]
                lambda _ev, did=did: self.domain_focus_requested.emit(int(did) + 1)
            )
            container.setCursor(Qt.PointingHandCursor)
            self._domains_layout.insertWidget(insert_at, container)
            self._domain_row_widgets.append(container)
            insert_at += 1

    def _build_kv_card(self, title: str, rows: list[tuple[str, str]]) -> QFrame:
        card, layout = make_card(title, margins=(18, 16, 18, 18), spacing=6)
        for label_text, key in rows:
            row = QHBoxLayout()
            row.setContentsMargins(0, 0, 0, 0)
            row.setSpacing(10)
            label = QLabel(label_text)
            label.setProperty("class", "label")
            label.setMinimumWidth(118)
            value = QLabel("—")
            value.setStyleSheet("font-weight: 700; color: #111827;")
            value.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
            value.setObjectName(f"inspector_value_{key}")
            row.addWidget(label)
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
