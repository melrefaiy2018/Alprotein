"""
pyqtgraph-based spectrum widget.

This is the Phase D1 deliverable: a near drop-in replacement for
:class:`SpectrumPlotWidget` that uses :mod:`pyqtgraph` instead of matplotlib.
pyqtgraph redraws at ~60 fps so we can drive live sliders without blitting
gymnastics. The matplotlib widget is kept for now so we can compare them
side-by-side before retiring it.

Public surface (intentionally a subset of the matplotlib widget's):

* signals: ``run_spectra_requested``, ``export_requested``, ``parameters_changed``
* methods: ``set_run_enabled``, ``get_parameters``, ``update_spectra``,
  ``update_spectrum_with_components``, ``clear``
"""

from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np
from PyQt5.QtCore import Qt, QTimer, pyqtSignal
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSizePolicy,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

try:
    import pyqtgraph as pg
    PYQTGRAPH_AVAILABLE = True
except ImportError:  # pragma: no cover — required by Phase D, defended for safety.
    pg = None  # type: ignore[assignment]
    PYQTGRAPH_AVAILABLE = False


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

ABS_PEN = (37, 99, 235)        # blue
FLU_PEN = (220, 38, 38)        # red
EXP_PEN = (107, 114, 128)      # gray (experimental)
STICK_PEN = (34, 197, 94, 160) # green, semi-transparent


def _configure_pg_defaults() -> None:
    """Set sensible pyqtgraph defaults once."""
    if not PYQTGRAPH_AVAILABLE:
        return
    pg.setConfigOption("background", "w")
    pg.setConfigOption("foreground", "k")
    pg.setConfigOption("antialias", True)


_configure_pg_defaults()


# ---------------------------------------------------------------------------
# Widget
# ---------------------------------------------------------------------------

class FastSpectrumWidget(QWidget):
    """pyqtgraph-backed spectrum panel with live controls."""

    export_requested = pyqtSignal()
    run_spectra_requested = pyqtSignal()
    parameters_changed = pyqtSignal(dict)

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self.wavelengths_abs: Optional[np.ndarray] = None
        self.absorption: Optional[np.ndarray] = None
        self.wavelengths_fl: Optional[np.ndarray] = None
        self.fluorescence: Optional[np.ndarray] = None
        self.stick_wavelengths: Optional[np.ndarray] = None
        self.stick_intensities: Optional[np.ndarray] = None
        self.experimental_x: Optional[np.ndarray] = None
        self.experimental_y: Optional[np.ndarray] = None

        # Debounce timer for live slider edits → ``parameters_changed``.
        # Slider drags fire ``valueChanged`` continuously; we coalesce them
        # so the workbench only triggers one recompute per drag.
        self._param_debounce = QTimer(self)
        self._param_debounce.setSingleShot(True)
        self._param_debounce.setInterval(350)
        self._param_debounce.timeout.connect(self._emit_parameters)

        self._setup_ui()

    # ------------------------------------------------------------------
    # UI
    # ------------------------------------------------------------------

    def _setup_ui(self) -> None:
        outer = QVBoxLayout(self)
        outer.setContentsMargins(8, 8, 8, 8)
        outer.setSpacing(8)

        header = QHBoxLayout()
        title = QLabel("Spectra (Fast)")
        title.setStyleSheet("font-size: 14px; font-weight: 600;")
        header.addWidget(title)
        header.addStretch(1)

        self._run_btn = QPushButton("▶ Run Spectra")
        self._run_btn.setProperty("class", "primary")
        self._run_btn.setEnabled(False)
        self._run_btn.clicked.connect(self.run_spectra_requested.emit)
        header.addWidget(self._run_btn)

        self._export_btn = QPushButton("📤 Export")
        self._export_btn.clicked.connect(self.export_requested.emit)
        header.addWidget(self._export_btn)
        outer.addLayout(header)

        # --- Controls row (live, debounced) ----------------------------
        controls = QFrame()
        controls.setProperty("class", "card")
        ctrl_grid = QGridLayout(controls)
        ctrl_grid.setContentsMargins(10, 8, 10, 8)
        ctrl_grid.setHorizontalSpacing(10)
        ctrl_grid.setVerticalSpacing(4)

        # Temperature
        self._temperature = QDoubleSpinBox()
        self._temperature.setRange(4.0, 400.0)
        self._temperature.setValue(300.0)
        self._temperature.setDecimals(0)
        self._temperature.setSuffix(" K")
        self._temperature.valueChanged.connect(self._schedule_param_emit)

        # Disorder σ
        self._disorder = QDoubleSpinBox()
        self._disorder.setRange(0.0, 800.0)
        self._disorder.setValue(130.0)
        self._disorder.setDecimals(0)
        self._disorder.setSuffix(" cm⁻¹")
        self._disorder.valueChanged.connect(self._schedule_param_emit)

        # Ensemble size
        self._ensemble = QSpinBox()
        self._ensemble.setRange(10, 5000)
        self._ensemble.setValue(300)
        self._ensemble.setSingleStep(10)
        self._ensemble.valueChanged.connect(self._schedule_param_emit)

        # x-axis unit
        self._x_axis = QComboBox()
        self._x_axis.addItems(["Wavelength (nm)", "Energy (cm⁻¹)", "Energy (eV)"])
        self._x_axis.currentIndexChanged.connect(self._refresh_plot)

        # Toggles
        self._show_abs = QCheckBox("Absorption")
        self._show_abs.setChecked(True)
        self._show_abs.stateChanged.connect(self._refresh_plot)

        self._show_flu = QCheckBox("Fluorescence")
        self._show_flu.setChecked(True)
        self._show_flu.stateChanged.connect(self._refresh_plot)

        self._show_sticks = QCheckBox("Sticks")
        self._show_sticks.setChecked(False)
        self._show_sticks.stateChanged.connect(self._refresh_plot)

        self._show_exp = QCheckBox("Experiment")
        self._show_exp.setChecked(False)
        self._show_exp.stateChanged.connect(self._refresh_plot)

        # CD / LD placeholders. The calculation isn't wired yet (P1.1
        # roadmap item); the toggle is here so the UI lands ahead of the
        # math and we can stub-render an empty curve.
        self._show_cd = QCheckBox("CD")
        self._show_cd.setEnabled(False)
        self._show_cd.setToolTip(
            "Circular dichroism — coming in roadmap item P1.1.\n"
            "Requires rotational strengths from the exciton eigenvectors."
        )
        self._show_ld = QCheckBox("LD")
        self._show_ld.setEnabled(False)
        self._show_ld.setToolTip(
            "Linear dichroism — coming in roadmap item P1.1.\n"
            "Requires a defined membrane normal."
        )

        # Layout: 2 rows × 4 columns of label/widget pairs
        def add(row: int, col: int, label: str, w: QWidget) -> None:
            ctrl_grid.addWidget(QLabel(label), row, col * 2)
            ctrl_grid.addWidget(w, row, col * 2 + 1)

        add(0, 0, "T", self._temperature)
        add(0, 1, "σ", self._disorder)
        add(0, 2, "N", self._ensemble)
        add(0, 3, "x", self._x_axis)
        ctrl_grid.addWidget(self._show_abs, 1, 0, 1, 1)
        ctrl_grid.addWidget(self._show_flu, 1, 1, 1, 1)
        ctrl_grid.addWidget(self._show_sticks, 1, 2, 1, 1)
        ctrl_grid.addWidget(self._show_exp, 1, 3, 1, 1)
        ctrl_grid.addWidget(self._show_cd, 1, 4, 1, 1)
        ctrl_grid.addWidget(self._show_ld, 1, 5, 1, 1)

        outer.addWidget(controls)

        # --- Experiment overlay card -----------------------------------
        exp_card = QFrame()
        exp_card.setProperty("class", "card")
        exp_layout = QHBoxLayout(exp_card)
        exp_layout.setContentsMargins(10, 6, 10, 6)
        exp_layout.setSpacing(10)

        load_btn = QPushButton("📂 Load experiment CSV…")
        load_btn.clicked.connect(self._on_load_experiment)
        exp_layout.addWidget(load_btn)

        self._exp_status = QLabel("No experiment loaded")
        self._exp_status.setProperty("class", "label")
        exp_layout.addWidget(self._exp_status, 1)

        self._rmse_label = QLabel("RMSE: —")
        self._rmse_label.setStyleSheet("font-family: monospace;")
        exp_layout.addWidget(self._rmse_label)

        clear_btn = QPushButton("✕")
        clear_btn.setFixedWidth(28)
        clear_btn.setToolTip("Clear experimental overlay")
        clear_btn.clicked.connect(self._on_clear_experiment)
        exp_layout.addWidget(clear_btn)

        outer.addWidget(exp_card)

        # --- Plot -------------------------------------------------------
        if PYQTGRAPH_AVAILABLE:
            self._plot = pg.PlotWidget()
            self._plot.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
            self._plot.showGrid(x=True, y=True, alpha=0.3)
            self._plot.setLabel("bottom", "Wavelength", units="nm")
            self._plot.setLabel("left", "Intensity (a.u.)")
            self._plot.addLegend(offset=(-10, 10))

            self._abs_curve = self._plot.plot(pen=pg.mkPen(color=ABS_PEN, width=2),
                                              name="Absorption")
            self._flu_curve = self._plot.plot(pen=pg.mkPen(color=FLU_PEN, width=2),
                                              name="Fluorescence")
            self._exp_curve = self._plot.plot(
                pen=pg.mkPen(color=EXP_PEN, width=1, style=Qt.DashLine),
                name="Experiment",
            )
            self._stick_item = pg.BarGraphItem(x=[], height=[], width=0.5,
                                               brush=pg.mkBrush(*STICK_PEN))
            self._plot.addItem(self._stick_item)
        else:
            self._plot = QLabel("pyqtgraph is not installed — install with `pip install pyqtgraph`.")
            self._plot.setAlignment(Qt.AlignCenter)
            self._plot.setStyleSheet("color: #8b2b2b; padding: 40px;")
            self._abs_curve = self._flu_curve = self._exp_curve = None
            self._stick_item = None

        outer.addWidget(self._plot, 1)

    # ------------------------------------------------------------------
    # Public API expected by the workbench
    # ------------------------------------------------------------------

    def set_run_enabled(self, enabled: bool) -> None:
        self._run_btn.setEnabled(enabled)

    def get_parameters(self) -> Dict[str, Any]:
        return {
            "temperature": float(self._temperature.value()),
            "disorder_sigma": float(self._disorder.value()),
            "n_ensemble": int(self._ensemble.value()),
            "include_vibronic": True,
            "use_ensemble": True,
        }

    def update_spectra(
        self,
        wavelengths_abs: np.ndarray,
        absorption: np.ndarray,
        wavelengths_fl: np.ndarray,
        fluorescence: np.ndarray,
    ) -> None:
        """Push freshly-computed absorption + fluorescence curves."""
        self.wavelengths_abs = np.asarray(wavelengths_abs)
        self.absorption = np.asarray(absorption)
        self.wavelengths_fl = np.asarray(wavelengths_fl)
        self.fluorescence = np.asarray(fluorescence)
        self._refresh_plot()
        self._update_rmse()

    def update_spectrum_with_components(
        self,
        wavelengths: np.ndarray,
        spectrum: np.ndarray,
        stick_wavelengths: Optional[np.ndarray] = None,
        stick_intensities: Optional[np.ndarray] = None,
        **_kwargs: Any,
    ) -> None:
        """Legacy compat path used by the simple-spectrum fallback."""
        self.wavelengths_abs = np.asarray(wavelengths)
        self.absorption = np.asarray(spectrum)
        if stick_wavelengths is not None and stick_intensities is not None:
            self.stick_wavelengths = np.asarray(stick_wavelengths)
            self.stick_intensities = np.asarray(stick_intensities)
        self._refresh_plot()

    def set_experimental(self, x: np.ndarray, y: np.ndarray, label: str = "") -> None:
        """Attach an experimental reference curve for overlay."""
        self.experimental_x = np.asarray(x, dtype=float)
        self.experimental_y = np.asarray(y, dtype=float)
        self._show_exp.setChecked(True)
        msg = label or f"{len(self.experimental_x)} points loaded"
        self._exp_status.setText(msg)
        self._refresh_plot()
        self._update_rmse()

    def clear_experimental(self) -> None:
        self.experimental_x = None
        self.experimental_y = None
        self._show_exp.setChecked(False)
        self._exp_status.setText("No experiment loaded")
        self._rmse_label.setText("RMSE: —")
        self._refresh_plot()

    def _on_load_experiment(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Load experimental spectrum",
            "",
            "CSV / text (*.csv *.dat *.txt);;All files (*)",
        )
        if not path:
            return
        try:
            data = np.loadtxt(path, delimiter=None, comments="#")
        except ValueError:
            # Try comma-delimited as a fallback.
            data = np.loadtxt(path, delimiter=",", comments="#")
        except OSError as e:
            self._exp_status.setText(f"Could not read: {e}")
            return
        if data.ndim != 2 or data.shape[1] < 2:
            self._exp_status.setText("Expected ≥ 2 numeric columns")
            return
        from pathlib import Path
        name = Path(path).name
        self.set_experimental(data[:, 0], data[:, 1], label=name)

    def _on_clear_experiment(self) -> None:
        self.clear_experimental()

    def _update_rmse(self) -> None:
        """Compute RMSE between calculated absorption and experiment (peak-normalized)."""
        if (
            self.experimental_x is None
            or self.experimental_y is None
            or self.wavelengths_abs is None
            or self.absorption is None
        ):
            self._rmse_label.setText("RMSE: —")
            return
        # Interpolate experiment onto the calculated wavelength grid (or vice
        # versa) and peak-normalize both before comparison.
        calc = self._normalize(self.absorption)
        exp_on_calc = np.interp(
            self.wavelengths_abs,
            self.experimental_x,
            self.experimental_y,
            left=np.nan,
            right=np.nan,
        )
        exp_on_calc = exp_on_calc / np.nanmax(np.abs(exp_on_calc)) if np.nanmax(np.abs(exp_on_calc)) > 0 else exp_on_calc
        mask = np.isfinite(exp_on_calc)
        if not np.any(mask):
            self._rmse_label.setText("RMSE: —")
            return
        rmse = float(np.sqrt(np.mean((calc[mask] - exp_on_calc[mask]) ** 2)))
        self._rmse_label.setText(f"RMSE: {rmse:.3f}")

    def clear(self) -> None:
        self.wavelengths_abs = None
        self.absorption = None
        self.wavelengths_fl = None
        self.fluorescence = None
        self.stick_wavelengths = None
        self.stick_intensities = None
        self.experimental_x = None
        self.experimental_y = None
        self._refresh_plot()

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _emit_parameters(self) -> None:
        self.parameters_changed.emit(self.get_parameters())

    def _schedule_param_emit(self, *_args) -> None:
        """Coalesce rapid slider edits — only emit after a quiet 350 ms."""
        self._param_debounce.start()

    def _x_transform(self, wl_nm: np.ndarray) -> np.ndarray:
        unit = self._x_axis.currentText()
        if unit.startswith("Wavelength"):
            return wl_nm
        # 10⁷ / λ_nm = cm⁻¹.
        with np.errstate(divide="ignore"):
            cm_inv = np.where(wl_nm > 0, 1.0e7 / wl_nm, 0.0)
        if "cm" in unit:
            return cm_inv
        # eV: hc / λ. λ in m → E in eV: 1239.84193 / λ_nm.
        return 1239.841984 / np.where(wl_nm > 0, wl_nm, np.nan)

    def _refresh_plot(self) -> None:
        if not PYQTGRAPH_AVAILABLE or self._abs_curve is None:
            return
        unit = self._x_axis.currentText()
        if unit.startswith("Wavelength"):
            self._plot.setLabel("bottom", "Wavelength", units="nm")
        elif "cm" in unit:
            self._plot.setLabel("bottom", "Energy", units="cm⁻¹")
        else:
            self._plot.setLabel("bottom", "Energy", units="eV")

        # Absorption
        if self._show_abs.isChecked() and self.wavelengths_abs is not None and self.absorption is not None:
            self._abs_curve.setData(
                self._x_transform(self.wavelengths_abs), self._normalize(self.absorption)
            )
        else:
            self._abs_curve.setData([], [])

        # Fluorescence
        if self._show_flu.isChecked() and self.wavelengths_fl is not None and self.fluorescence is not None:
            self._flu_curve.setData(
                self._x_transform(self.wavelengths_fl), self._normalize(self.fluorescence)
            )
        else:
            self._flu_curve.setData([], [])

        # Sticks
        if (
            self._show_sticks.isChecked()
            and self.stick_wavelengths is not None
            and self.stick_intensities is not None
        ):
            xs = self._x_transform(self.stick_wavelengths)
            heights = self._normalize(self.stick_intensities)
            self._stick_item.setOpts(x=xs, height=heights, width=0.4)
        else:
            self._stick_item.setOpts(x=[], height=[], width=0.4)

        # Experiment
        if self._show_exp.isChecked() and self.experimental_x is not None and self.experimental_y is not None:
            self._exp_curve.setData(
                self._x_transform(self.experimental_x), self._normalize(self.experimental_y)
            )
        else:
            self._exp_curve.setData([], [])

    @staticmethod
    def _normalize(y: np.ndarray) -> np.ndarray:
        y = np.asarray(y, dtype=float)
        peak = float(np.nanmax(np.abs(y))) if y.size else 0.0
        if peak <= 0.0 or not np.isfinite(peak):
            return y
        return y / peak


__all__ = ["FastSpectrumWidget", "PYQTGRAPH_AVAILABLE"]
