# Alprotein — Audit & Improvement Roadmap

**Version reviewed:** 0.2.0 · **Date:** 2026-05-14 · **Scope:** capabilities, features, UI/UX

This roadmap is grounded in a direct read of the source tree (~21k LOC) on `main`. Items are prioritized P0 (do first) → P3 (nice-to-have). Effort is rough: S = ≤1 day, M = 2–5 days, L = 1–2 weeks, XL = >2 weeks.

---

## 1. Executive Summary

Alprotein is a Python + PyQt5 toolkit for **excitonic Hamiltonians and steady-state optical spectra of pigment–protein complexes** (chlorophylls, BChl, pheophytin). The scientific core is solid: CDC site energies, TrEsp / point-dipole couplings, ensemble-averaged absorption & fluorescence with vibronic coupling under a double-overdamped Brownian-oscillator spectral density (Renger-style), with on-disk caching of g(t).

The biggest immediate liabilities are **not scientific** — they are **product hygiene**: three competing GUI stacks coexist in `Alprotein/gui/`, a research `scratch/` directory ships inside the package, four `.py` modules are empty stubs, and the public `launch_gui()` entry point still points at the legacy window instead of the Scientific Workbench that `tools/launch_workbench.py` launches. Fix those before adding features.

The biggest scientific opportunities are **CD/LD spectra** (nearly free given existing transition dipoles), **dynamics** (Redfield / modified Redfield → 2D electronic spectra), and **MD-trajectory disorder** sampling. These are what would move the package from "static-spectra calculator" to "publication-quality exciton suite."

---

## 2. What's Here Today (Code Map)

| Area | Files | Notes |
|---|---|---|
| Core | `core/protein_structure.py` (906), `core/pigment_system.py` (570), `core/abstract_pigments.py` (837), `core/constants.py`, `core/enhanced_atom.py` | Clean, well-scoped. |
| Calculators | `calculators/site_energy_calculator.py` (826), `coupling_calculator.py` (325), `hamiltonian_calculator.py` (219), `exciton_calculator.py` (320), `spectra_calculator.py` (944) | Spectra calc has disk-cached g(t) — good. |
| GUI (live) | `gui/workbench_app.py` (520) + `gui/workbench_window.py` (1796) + ~17 widgets | Workbench is the real product. `workbench_window.py` is a god-object. |
| GUI (legacy) | `gui/app.py`, `app_pro.py`, `main_window.py`, `main_window_pro.py` | Dead weight. `gui/__init__.launch_gui` still imports the legacy `AlproteinGUI`. |
| 3D viewer | `gui/widgets/protein_viewer.py` via **py3Dmol + QWebEngineView** | Already modern; not the custom-matplotlib I expected. |
| Utils | `utils/cdc_analysis.py` (698), `utils/calculate_exciton_distribution.py` (681), `utils/atom_level_cdc.py`, `utils/pdb_writer.py`, `utils/calculations.py` | Active. |
| Empty stubs | `utils/export_utils.py`, `utils/gui_utils.py`, `visualization/energy_plots.py`, `visualization/spatial_analysis.py` | 0 bytes — delete or implement. |
| Scratch (shipped in source tree) | `Alprotein/scratch/WSCP/Q57D/...` | Excluded from `packages.find` but lives inside the package dir. Move out. |
| Tests | 19 `test_*.py` files | Many are smoke/ad-hoc scripts, not proper pytest. |
| Docs | `docs/CHANGELOG.md`, `docs/PR.md`, `docs/gui-fixes/*.md` | Internal fix notes, not user docs. No Sphinx build. |

---

## 3. Bugs / Footguns Spotted

| # | Severity | Issue | Fix |
|---|---|---|---|
| B1 | **HIGH** | `Alprotein/gui/__init__.py` exports `launch_gui()` that constructs `AlproteinGUI` (legacy), while `tools/launch_workbench.py` launches the Scientific Workbench. Two different "official" UIs depending on entry point. | Point `launch_gui` at `workbench_app.main`. |
| B2 | MED | Four empty `.py` files (`export_utils.py`, `gui_utils.py`, `energy_plots.py`, `spatial_analysis.py`) are importable and confusing. | Delete or implement. |
| B3 | MED | `Alprotein/scratch/` is inside the installed package directory. `pyproject.toml` excludes it from wheel, but it pollutes the source tree and IDE search. | Move to a top-level `scratch/` outside `Alprotein/`. |
| B4 | MED | `README.md` author emails are placeholders (`@example.edu`) — same in `pyproject.toml`. | Fix before any release. |
| B5 | MED | `pyproject.toml` URLs point to `Alprotein-Alpha` but the repo is `AlProtein`. Broken links. | Update URLs. |
| B6 | LOW | Conflicting GUI installation: `gui/__init__.py` does `from .app import AlproteinGUI` at import time. If PyQt is missing, the **package-level** `try/except ImportError` in `Alprotein/__init__.py` will hide every other GUI import error too, not just missing PyQt. | Narrow the `try` to PyQt only, or lazy-import. |
| B7 | LOW | `Alprotein/__init__.py` `quick_calculation()` calls `calculator.construct_hamiltonian(kwargs)` — `kwargs` is a dict, but `construct_hamiltonian` expects `params=...`. Inspect signature. | Verify and pass `params=kwargs`. |
| B8 | LOW | Default Python is 3.8 in classifiers but `conda create -n alprotein python=3.11` in README. NumPy ≥1.18 / scipy ≥1.6 are very old floors. | Bump min to 3.10. |
| B9 | LOW | `_print_welcome()` exists but is commented-out — dead code. | Remove. |

---

## 4. Priority 0 — Stabilize Before Anything New (≈1 week)

**Goal:** one canonical UI, clean source tree, working public entry points, a real test suite.

### P0.1 — Collapse the three GUI stacks (M)
- Keep `workbench_*` and `gui/widgets/*` and `gui/dialogs/*`.
- Delete `app.py`, `app_pro.py`, `main_window.py`, `main_window_pro.py` and any widgets only they reference.
- Fix `gui/__init__.py:launch_gui` to call `workbench_app.main`.
- Add `python -m Alprotein` entry point (`__main__.py`) that launches the workbench.
- Register a console script in `pyproject.toml`: `alprotein = Alprotein.gui.workbench_app:main`.

### P0.2 — Clean source tree (S)
- Remove empty stubs (B2) or implement them in P1/P2.
- Move `Alprotein/scratch/` → top-level `scratch/`.
- Fix author emails (B4) and repo URLs (B5).
- Remove dead `_print_welcome` (B9), `quick_calculation` bug (B7).

### P0.3 — Refactor the 1796-line workbench window (M)
- `workbench_window.py` mixes UI construction, signal wiring, calculation orchestration, plot picking, and export. Extract into:
  - `controllers/calculation_controller.py` — `CalculationWorker` + dispatch
  - `controllers/results_controller.py` — `on_calculation_finished`, plotting glue
  - `controllers/export_controller.py` — Hamiltonian/spectrum/table exports
  - `workbench_window.py` becomes a thin assembler (<400 lines)
- Same treatment for `spectrum_widget.py` (939 lines) if it has comparable scope creep.

### P0.4 — Real pytest suite (M)
- Audit the 19 `test_*.py`. Many are import smoke tests or GUI launch scripts. Keep them if useful, but add:
  - `tests/unit/` — pure-Python tests for `CouplingCalculator`, `HamiltonianCalculator`, `SpectraCalculator` against known toy systems (dimer with analytic eigenvalues, monomer spectrum check).
  - `tests/regression/` — golden-file comparison for CP29, IsiA, CP24 (already partially present) with `pytest.approx` and JSON-stored references.
  - GUI tests stay isolated under `tests/gui/` and are skipped when `QT_QPA_PLATFORM=offscreen` is unset on CI.
- Set up `pytest-cov`, target ≥70% on non-GUI code.

### P0.5 — CI (S)
- GitHub Actions: lint (`ruff` or keep `flake8`/`black`/`isort`), `pytest -q`, build wheel, docs build. Test matrix: 3.10 / 3.11 / 3.12 on ubuntu + macos.

---

## 5. Priority 1 — High-Leverage Scientific Capabilities (≈3–4 weeks)

These give the biggest scientific return per unit effort.

### P1.1 — Circular Dichroism (CD) & Linear Dichroism (LD) spectra (S)
**Why:** You already have exciton eigenvectors, site positions, and transition dipoles. Adding CD is ~50 lines of code and is the standard companion to absorption for chiral pigment arrays. LD is similar with a defined membrane normal.
- New method `SpectraCalculator.calculate_cd_spectrum(...)` using the rotational strength
  R_k = (π / λ_k) · Σ_{i<j} (R_i − R_j) · (μ_i × μ_j) · c_{ik} c_{jk}.
- Add toggle in spectrum widget; reuse the existing line-broadening machinery.

### P1.2 — Site-energy fitting / inverse mode (M)
**Why:** This is the #1 thing experimental groups want — fit site energies to a measured absorption (and optionally CD) curve while holding couplings/structure fixed.
- Wrap `SpectraCalculator.calculate_spectrum` with `scipy.optimize.least_squares` or `differential_evolution`.
- GUI tab: load experimental spectrum (CSV / two-column), select which site energies are free, bounds, run with live plot.
- Output: best-fit energies + parameter covariance.

### P1.3 — Charge-transfer (CT) state support (M)
**Why:** Reaction-center modelers cannot use Alprotein today because CT states aren't representable.
- Extend Hamiltonian builder to accept user-defined CT basis states with energies and mixing matrix elements to neighboring locally-excited states.
- Add CT entries to the Hamiltonian table widget.

### P1.4 — Modified-Redfield dynamics & 2D electronic spectra (L)
**Why:** This is the obvious next scientific frontier given you already compute the same g(t) used in Redfield rates.
- Add `calculators/redfield_calculator.py` for population transfer rates between exciton states.
- Add `calculators/two_d_calculator.py` for rephasing/non-rephasing 2DES with the response-function formalism.
- New "Dynamics" tab in the workbench with population traces and 2D plots.
- Caveat: this is publication-grade work; budget time accordingly.

### P1.5 — MD-trajectory ensemble (M)
**Why:** Replaces ad-hoc Gaussian disorder (σ in cm⁻¹) with **structurally sampled** disorder, which is the modern standard.
- `ProteinStructure.from_trajectory(top, traj, stride=...)` using MDAnalysis (optional dep).
- Compute site energies & couplings per frame, build an ensemble of Hamiltonians, average the spectra.
- Show "MD ensemble" as an alternative to "Gaussian disorder" in the spectra widget.

### P1.6 — Transition Density Cube (TDC) coupling option (M)
**Why:** TDC is the reference standard for short-range couplings where TrEsp under-converges.
- Read precomputed cube files (`.cube`) for each pigment type.
- Use Coulomb integration on the grid.
- Make it a third method in the existing TrEsp / point-dipole selector.

### P1.7 — Pigment library expansion (S)
**Why:** Currently only ChlA / ChlB / Pheophytin. Adding the rest unlocks more systems.
- BChl a, BChl b, Chl c, Chl d, Chl f, carotenoids (S1/S2 with caveats), phycobilins.
- Provide TrEsp / CDC parameter sets from the literature with citations in docstrings.

---

## 6. Priority 2 — UI/UX (≈2–3 weeks, in parallel with P1)

The workbench is functional but feels like a research tool, not a product. Treat the user as a biophysicist who wants to load a PDB and *see things*.

### P2.1 — Project save / load (M)
- `.alproj` file: PDB path + parameters + computed Hamiltonian + spectra (HDF5 inside a zip).
- Recent files menu.
- Eliminates "redo everything when I reopen."

### P2.2 — Live parameter sliders (M)
- Sliders for T, σ-disorder, dielectric → recompute spectrum on release (or with debounce).
- Use **pyqtgraph** for the spectrum widget instead of matplotlib (5–10× faster redraw, native zoom/pan).
- Keep matplotlib only for "publication export" mode.

### P2.3 — Pigment ↔ Hamiltonian ↔ Spectrum cross-highlighting (M)
- Click pigment in 3D viewer → highlight row/col in Hamiltonian + contributing exciton band in spectrum.
- Click an exciton band → highlight pigments weighted by |c_ik|² on the 3D structure (color heatmap).
- Click a Hamiltonian off-diagonal → draw a line between the two pigments in 3D.
- This is the killer demo feature.

### P2.4 — Domain visualization (S)
- The build_domains step already groups pigments by coupling. Color them by domain in the 3D viewer.

### P2.5 — Progress + cancellation (S)
- Long calculations (ensemble spectra) need a real progress bar with a cancel button. The `CalculationWorker` is a QThread — wire a stop flag through the ensemble loop.

### P2.6 — In-app docs / tooltips with references (S)
- Each parameter (E_0a, σ-disorder, f_val, dielectric_cdc, …) needs a tooltip with units, typical range, and a citation. Right now users have to guess.

### P2.7 — Dark mode + theme polish (S)
- Centralize colors in `gui/styles.py`. Add a toggle.

### P2.8 — Onboarding (S)
- First-run dialog: "Try with example PDB (CP29 / IsiA / WSCP)" buttons that load bundled files.

### P2.9 — Export polish (S)
- "Publication figure" export: matplotlib-rendered PNG/SVG at user-chosen DPI, with proper axes, legend, and a single-call style.
- HDF5 export for Hamiltonian + spectra together.

---

## 7. Priority 3 — Nice to Have

- **Web/Notebook embedding.** A `Alprotein.notebook` thin API for Jupyter (already feasible — just document the patterns).
- **Flask/FastAPI backend** (the `api` optional dep already lists flask but no code exists). Useful for headless HPC use.
- **GPU-accelerated ensemble averaging** via `cupy` or `jax` for the disorder loop in `SpectraCalculator`.
- **Plugin system** for custom pigment classes and custom couplings.
- **CLI** (`alprotein run config.yaml`) for HPC batch jobs.
- **Sphinx docs site** on Read the Docs.
- **Hosted example gallery** with notebooks for CP29, FMO, LH2, PSII core, IsiA.

---

## 8. Suggested 90-Day Plan

| Weeks | Focus | Deliverables |
|---|---|---|
| 1 | P0 stabilize | One GUI, clean tree, fixed entry points, CI green |
| 2 | P0 refactor + tests | Split workbench_window; pytest ≥70% non-GUI |
| 3 | P1.1 + P1.7 | CD spectra shipped; 3 new pigment classes |
| 4–5 | P1.2 | Site-energy fitter (tab + scipy + tests) |
| 6 | P2.1 + P2.5 + P2.6 | Project save/load, progress/cancel, tooltips |
| 7 | P2.2 + P2.9 | pyqtgraph spectrum + publication export |
| 8 | P2.3 | Cross-highlighting (pigment ↔ Hamiltonian ↔ spectrum) |
| 9–11 | P1.5 + P1.6 | MD ensembles + TDC couplings |
| 12 | Release 0.3.0 | Tag, docs, PyPI, conda-forge submission |

Then 0.4.0 = P1.3 (CT) + P1.4 (Redfield/2DES) over the next quarter.

---

## 9. What I'd Do First If I Had One Day

1. Fix `gui/__init__.py:launch_gui` to point at the workbench (10 min).
2. Add `[project.scripts] alprotein = Alprotein.gui.workbench_app:main` (5 min).
3. Delete the four empty stub files (2 min).
4. Move `Alprotein/scratch/` out of the package (10 min).
5. Implement P1.1 (CD spectra) — single afternoon, immediate scientific value.
6. Add the bundled-example dropdown in the file loader (P2.8) — gives reviewers a 30-second "wow" path.

That alone takes the project from "research code with a GUI bolted on" to "tool a colleague could install and use in an afternoon."
