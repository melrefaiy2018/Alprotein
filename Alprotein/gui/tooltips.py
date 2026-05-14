"""
Central registry of parameter tooltips with units, ranges, and citations.

Use :func:`apply_tooltips` to attach rich tooltips to a mapping of
``parameter_key -> QWidget``. The tooltip text is rendered as Qt rich text so
it works inside the standard ``QToolTip``.

Each entry can be referenced from a help panel or the command palette to give
users a "What is this?" lookup without leaving the workbench.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Mapping, Optional

from PyQt5.QtWidgets import QWidget


@dataclass(frozen=True)
class ParameterDoc:
    """Documentation entry for a single calculation parameter."""

    label: str
    description: str
    units: str = ""
    typical_range: str = ""
    citation: str = ""
    doi: str = ""

    def as_html(self) -> str:
        """Render the entry as an HTML tooltip."""
        parts = [f"<b>{self.label}</b>", f"<p>{self.description}</p>"]
        meta = []
        if self.units:
            meta.append(f"<b>Units:</b> {self.units}")
        if self.typical_range:
            meta.append(f"<b>Typical:</b> {self.typical_range}")
        if meta:
            parts.append("<br>".join(meta))
        if self.citation:
            cite = self.citation
            if self.doi:
                cite = f"{cite} · <code>{self.doi}</code>"
            parts.append(f'<p style="color:#6b7280; font-size: 11px;">{cite}</p>')
        return "".join(parts)


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

PARAMETER_DOCS: Dict[str, ParameterDoc] = {
    "dielectric_cdc": ParameterDoc(
        label="ε(CDC)",
        description=(
            "Effective dielectric constant used in the Coulomb sum that produces "
            "site-energy shifts via the charge-density coupling (CDC) method. "
            "Higher ε screens electrostatic shifts more strongly."
        ),
        units="dimensionless",
        typical_range="1.5 – 4 (interior pigments)",
        citation="Adolphs & Renger, Biophys. J. 91, 2778 (2006)",
        doi="10.1529/biophysj.105.079483",
    ),
    "dielectric_tresp": ParameterDoc(
        label="ε(TrEsp)",
        description=(
            "Dielectric screening applied to TrEsp couplings between transition "
            "charges of different pigments. Usually left at 1 because environment "
            "screening is absorbed into the oscillator-strength factor f_osc."
        ),
        units="dimensionless",
        typical_range="1.0 – 2.0",
        citation="Madjet, Abdurahman & Renger, J. Phys. Chem. B 110, 17268 (2006)",
        doi="10.1021/jp0615398",
    ),
    "f_osc": ParameterDoc(
        label="Oscillator strength factor (f_osc)",
        description=(
            "Empirical scaling factor applied to TrEsp transition charges to "
            "match the in-protein oscillator strength of the Qy transition."
        ),
        units="dimensionless",
        typical_range="0.6 – 0.9 (Chl a Qy)",
        citation="Madjet et al., J. Phys. Chem. B 110, 17268 (2006)",
        doi="10.1021/jp0615398",
    ),
    "e0a": ParameterDoc(
        label="Reference site energy — Chl a (E₀,a)",
        description=(
            "Reference Qy transition energy of monomeric chlorophyll a in vacuum, "
            "before adding the CDC electrostatic shift. CDC shifts are evaluated "
            "relative to this value."
        ),
        units="cm⁻¹",
        typical_range="14 700 – 15 100",
        citation="Renger & Schlodder, J. Photochem. Photobiol. B 104, 126 (2011)",
        doi="10.1016/j.jphotobiol.2011.02.016",
    ),
    "e0b": ParameterDoc(
        label="Reference site energy — Chl b (E₀,b)",
        description=(
            "Reference Qy transition energy of monomeric chlorophyll b in vacuum."
        ),
        units="cm⁻¹",
        typical_range="15 400 – 15 900",
        citation="Linnanto & Korppi-Tommola, Phys. Chem. Chem. Phys. 8, 663 (2006)",
        doi="10.1039/B513086G",
    ),
    "temperature": ParameterDoc(
        label="Temperature (T)",
        description=(
            "Sample temperature used in the Bose–Einstein occupation of phonon "
            "modes when evaluating the lineshape function g(t). Drives the "
            "homogeneous linewidth and the Stokes shift."
        ),
        units="K",
        typical_range="77 K (cryogenic) – 300 K (room)",
    ),
    "disorder_sigma": ParameterDoc(
        label="Static disorder σ",
        description=(
            "Standard deviation of the Gaussian distribution from which site-"
            "energy realisations are drawn for the ensemble average. Larger σ "
            "produces broader, more featureless spectra."
        ),
        units="cm⁻¹ (FWHM in many publications — confirm convention)",
        typical_range="80 – 200",
    ),
    "n_ensemble": ParameterDoc(
        label="Ensemble size (N)",
        description=(
            "Number of disorder realisations averaged when computing absorption "
            "/ fluorescence. Convergence is typically reached around 200–500."
        ),
        units="count",
        typical_range="200 – 1000",
    ),
    "coupling_method": ParameterDoc(
        label="Coupling method",
        description=(
            "How pigment–pigment couplings are evaluated. <b>TrEsp</b> uses "
            "transition charges fitted to the QM transition density. "
            "<b>Point dipole</b> is the leading multipole approximation; it "
            "breaks down at short separations (≲10 Å)."
        ),
        citation="Renger, Photosynth. Res. 102, 471 (2009)",
        doi="10.1007/s11120-009-9472-9",
    ),
    "domain_cutoff": ParameterDoc(
        label="Domain coupling cutoff",
        description=(
            "Couplings stronger than this threshold define edges of the exciton "
            "domain graph. Domains are diagonalised independently when building "
            "the spectrum, which is much faster than diagonalising the full "
            "Hamiltonian for large systems."
        ),
        units="cm⁻¹",
        typical_range="15 – 30",
    ),
}


# ---------------------------------------------------------------------------
# Application helpers
# ---------------------------------------------------------------------------

def get_doc(key: str) -> Optional[ParameterDoc]:
    """Return the documentation entry for ``key`` or ``None``."""
    return PARAMETER_DOCS.get(key)


def tooltip_html(key: str) -> str:
    """Return the HTML tooltip for ``key`` (empty string if unknown)."""
    doc = get_doc(key)
    return doc.as_html() if doc else ""


def apply_tooltips(mapping: Mapping[str, QWidget]) -> None:
    """Attach rich tooltips to widgets keyed by parameter name."""
    for key, widget in mapping.items():
        html = tooltip_html(key)
        if html and widget is not None:
            widget.setToolTip(html)


def known_keys() -> Iterable[str]:
    """Return the registered parameter keys."""
    return PARAMETER_DOCS.keys()


__all__ = [
    "ParameterDoc",
    "PARAMETER_DOCS",
    "get_doc",
    "tooltip_html",
    "apply_tooltips",
    "known_keys",
]
