# UI/UX Audit — Mockup vs. Shipped Code

**Mockup:** `docs/ux-mockups/workbench_v3.html`
**Branch:** `feat/gui-ux-overhaul` · 8 commits ahead of `main`
**Date:** 2026-05-14

Status key:
- ✅ **Done** — visible and functional in the app right now
- 🟡 **Partial** — wired in code but not visible / not finished
- ❌ **Missing** — nothing in the code corresponds

I went pane-by-pane through the mockup. Each row is one visible thing in the mockup.

---

## Section 1 — Full Window, Default Layout

### 1a. Title bar

| Mockup element | Status | Where it actually is | Gap |
|---|---|---|---|
| Window title shows project name + dirty marker (`"CP29_Doran.alproj • modified"`) | 🟡 | Title is set on save/open in `workbench_app._save_to` / `load_project_path`. Dirty marker tracked in `_project_dirty` but **never rendered**. | Add `• modified` suffix in the title when `_project_dirty=True`. |
| Subtitle muted text (`v0.3.0 · py3Dmol · pyqtgraph`) | ❌ | Nothing. | Cosmetic only — easy add. |

### 1b. Menu bar

| Mockup element | Status | Where it actually is | Gap |
|---|---|---|---|
| File, Edit, Structure, Calculate, Analyze, Visualize, Window, Help | 🟡 | All present except **Window**. Visualize is now mostly tab-navigation. | "Window" menu would hold show/hide Inspector etc. — not implemented. |

### 1c. Toolbar row (under the menu)

| Mockup element | Status | Where it actually is | Gap |
|---|---|---|---|
| 📂 Open / 💾 Save / 🕒 Recent ▾ | 🟡 | Lives in File menu and ⌘K. | **No visible toolbar.** This is the single biggest reason you said it looks different. |
| ↶ Undo / ↷ Redo | 🟡 | Edit menu + ⌘Z. | No visible toolbar. |
| ▶ Run All / ■ Stop | 🟡 | Calculate menu + ⌘. | No visible toolbar. |
| 📤 Export ▾ / 🎨 Theme: Light ▾ | 🟡 | File menu + View menu. | No visible toolbar. |
| ⌘K Command palette button at the right | 🟡 | Bound to ⌘K shortcut + View menu item. | No visible button. |
| Autosave chip ("Auto-save ✓") | ❌ | No auto-save timer in code. | Not implemented at all. |

### 1d. Left sidebar — Project card

| Mockup element | Status | Where it actually is | Gap |
|---|---|---|---|
| Drag-drop "🧬 Drag .pdb or .alproj here" | 🟡 | `DragDropArea` in `tools_panel.py` accepts files, but only `.pdb` (the workbench routes `.alproj` through `File → Open Project…` instead). | Make drop area accept `.alproj` too. |
| Project name + Pigments / Atoms / Validated stats | ✅ | `tools_panel.update_project_info`. | — |
| Auto-save badge | ❌ | No auto-save. | — |

### 1e. Left sidebar — Settings card (mockup shows **sliders**)

| Mockup element | Status | Where it actually is | Gap |
|---|---|---|---|
| ε(CDC), ε(TrEsp), f_osc, E₀(CLA), E₀(CLB) | ✅ | `tools_panel.create_settings_section`. | Rendered as **spin boxes**, not sliders. Functionally equivalent but visually different. |
| Coupling method dropdown (TrEsp / Point dipole / TDC planned) | 🟡 | Not in `tools_panel` — lives inside `hamiltonian_widget`. | Cosmetic placement difference. |
| T, σ disorder, N ensemble sliders | 🟡 | Live in `hamiltonian_widget` and the new `FastSpectrumWidget`, not in the left settings card. | Mockup puts everything in one card; we have them split per-tab. |
| "Apply & Recalculate" primary button | ❌ | No global Apply — every spin box is live. | Different interaction model, intentional. |
| Tooltip with citation on hover | ✅ | `tooltips.py` registry, `apply_tooltips()` in `tools_panel`. | — |

### 1f. Left sidebar — Results card

| Mockup element | Status | Where it actually is | Gap |
|---|---|---|---|
| Site energies row with status dot + summary text | ✅ | `tools_panel.update_result_status`. | — |
| Hamiltonian row | ✅ | Same. | — |
| Spectrum row | ✅ | Same. | — |
| Exciton distribution row | ✅ | Same. | — |
| "View All Results" button | 🟡 | `create_results_section` defines a `view_all_btn` but the **section is never added to the sidebar layout** (it's defined-and-orphaned in `tools_panel.py`). | Either wire it in or delete it. |

### 1g. Centre — Tab bar

| Mockup element | Status | Where it actually is | Gap |
|---|---|---|---|
| 🔬 3D Structure | ✅ | Tab 0. | — |
| 📐 Hamiltonian | ✅ | Tab 1. | — |
| 📊 Spectra | ✅ | Tab 2 (matplotlib, legacy). | — |
| 📋 Data & Analysis | ✅ | Tab 4. | — |
| 🎯 Fit to Experiment **NEW** | ❌ | Phase E — not started. | — |
| (extra) ⚡ Spectra (Fast) | n/a | Tab 3, added for A/B comparison with the matplotlib widget. | Will be removed once we delete the legacy spectrum widget. |

### 1h. Centre — Hamiltonian heatmap mock

| Mockup element | Status | Where it actually is | Gap |
|---|---|---|---|
| Title "Excitonic Hamiltonian · 13 × 13 · units: cm⁻¹" | 🟡 | Different layout — `hamiltonian_widget` shows a "site energy shifts" bar chart + a table, not a heatmap. | **No heatmap exists.** Mockup heatmap is genuinely missing. |
| Symmetrize toggle, Log scale, Copy, Export buttons | 🟡 | Export exists, others don't. | — |
| 13×13 colored grid | ❌ | — | — |
| Selected-Pair card (J, distance, angle, method) | ❌ | — | — |
| Eigenstates list | 🟡 | Eigenvalues are computed and stored on `workbench.eigenvalues`. The "View Full Matrix" button shows them in a dialog. No always-visible list. | — |

### 1i. Right sidebar — Inspector

| Mockup element | Status | Where it actually is | Gap |
|---|---|---|---|
| Selection title (e.g. "CLA-602") | ✅ | `inspector_panel.set_pigment`. | — |
| Chain/resi, Site E, ΔE shift, μ, Domain, Top k, Top partner | ✅ | All wired through `snapshot_for_pigment`. | μ shows "—" until pigments expose a dipole attribute. |
| Pin / Lock energy buttons | ❌ | No pin/lock concept. | — |
| Override E spin box | ✅ | Apply / Reset buttons. | — |

### 1j. Right sidebar — Live View Tools card

| Mockup element | Status | Where it actually is | Gap |
|---|---|---|---|
| Color-by dropdown (Pigment type / Site energy / \|c\|²(k) / Domain) | 🟡 | The exciton scrub bar paints by \|c\|²; domain colors are applied separately when domains update. **No single dropdown to pick one of the four modes.** | — |
| "Show couplings ≥" slider | ❌ | — | — |
| μ vectors / Mg-N₄ plane / Labels checkboxes | ❌ | — | — |

### 1k. Right sidebar — Notebook card

| Mockup element | Status | Where it actually is | Gap |
|---|---|---|---|
| Free-text notes (saved in .alproj) | 🟡 | `Project.notebook` field exists and is round-tripped. **No widget** to edit it. | Add a QPlainTextEdit. |

### 1l. Status bar

| Mockup element | Status | Where it actually is | Gap |
|---|---|---|---|
| Running message ("Computing spectrum — ensemble 142 / 300") | ✅ | `update_status` shows progress messages. | — |
| Progress bar | ✅ | `setup_status_bar`. | — |
| Cancel button | ✅ | Added in Phase A. | Only visible while a calc is running. |
| Elapsed / ETA | 🟡 | Workers don't report ETA today. | — |
| Cache hit indicator | ❌ | g(t) cache is silent. | — |
| RAM | ❌ | — | — |

### 1m. Error toast (bottom-right)

| Mockup element | Status | Where it actually is | Gap |
|---|---|---|---|
| Severity-coloured toast with "Copy details ›" | ✅ | `widgets/toast.py` + `ToastManager`; wired to `calculation_error`. | — |

---

## Section 2 — 3D Structure Tab

| Mockup element | Status | Where it actually is | Gap |
|---|---|---|---|
| py3Dmol scene | ✅ | `Viewer3D` (already used py3Dmol pre-this-branch). | — |
| Pigment click → callout + Python notification | ✅ | C2: QWebChannel bridge → `pigment_selected`. The mockup shows a coloured callout "CLA-602" floating over the selected pigment — we get the 3Dmol label (the existing one) plus Inspector update. | Mockup-style floating callout in the workbench's own theme not built. |
| 🎯 Center / 📷 Snapshot / ⟳ Reset overlay buttons | 🟡 | "Reset 3D View" exists in Visualize menu + Ctrl+R. No on-canvas overlay row. | — |
| Bottom overlay caption ("Color: \|c\|² of exciton k=3") | 🟡 | The readout next to the scrub slider shows `k = N / N · E = … cm⁻¹` but not "Color: …". | — |
| Exciton k scrub slider | ✅ | C3: `ProteinViewer._scrub_slider`. | — |
| Right sidebar — Inspector | ✅ | Same Inspector panel as Section 1i. | — |
| Right sidebar — Domains card with coloured swatches and a per-domain pigment count | 🟡 | Domains are tracked in `workbench._domains`; Inspector shows "Domain Dx" for the selected pigment. **No dedicated Domains card.** | — |
| Coupling cutoff slider in Domains card | 🟡 | Lives in `hamiltonian_widget`, not in the Domains card. | Placement difference. |

---

## Section 3 — Spectra Tab (Live sliders + experiment overlay)

| Mockup element | Status | Where it actually is | Gap |
|---|---|---|---|
| Spectrum plot (pyqtgraph) | ✅ | `FastSpectrumWidget` in tab "⚡ Spectra (Fast)". | The mockup shows it as *the* Spectra tab; we have it as an extra tab while the old matplotlib one still exists. |
| Legend (Absorption, Fluorescence, Experiment, Stick eigenvalues) | ✅ | pyqtgraph `addLegend`. | — |
| Per-exciton oscillator strength bar plot below the main plot | ❌ | Not in `FastSpectrumWidget`. | — |
| Right sidebar — Live Controls card (T, σ, N, vibronic, sticks, CD, LD, x-axis) | ✅ | Same controls inline at the top of the Fast tab. CD/LD are disabled placeholders. | Mockup puts them in a right-sidebar card; we put them above the plot. Functionally equivalent. |
| Right sidebar — Experiment Overlay card (Load CSV, Scale ×, Shift, RMSE, "Fit site energies →") | 🟡 | Load CSV ✅, RMSE ✅; **no Scale or Shift inputs; no "Fit" button** (that's Phase E). | — |

---

## Section 4 — Cross-cutting Patterns

| Mockup element | Status | Where it actually is | Gap |
|---|---|---|---|
| Project model (.alproj) | ✅ | `Alprotein/gui/project.py`. | — |
| Recent files | ✅ | `Alprotein/gui/recent_files.py`. | — |
| Auto-save every 60 s | ❌ | No timer. | — |
| Undo / Redo with visible "history panel" | 🟡 | `Alprotein/gui/undo.py`, Edit menu, palette. **No visible history list.** | — |
| Progress + cancellation | ✅ | Status bar Cancel button + `CalculationCancelled`. | — |
| In-app docs / tooltips with citations | 🟡 | Tooltips with DOIs ✅. **No side panel** rendering full markdown docs ("Help → What is…?"). | — |
| Command palette ⌘K | ✅ | `widgets/command_palette.py`; 25 commands registered. | — |
| Theme system (Light / Dark) | ✅ | `Alprotein/gui/theme.py` · View → Theme · Ctrl+Shift+T. | — |
| Error toasts | ✅ | `widgets/toast.py`. | — |
| Keyboard shortcuts | ✅ | Wired in menu actions + palette. | Some mockup-only shortcuts (e.g. ⌘. as Cancel) are present. |

---

## Section 5 — What This Replaces / Removes

| Mockup intent | Status | Notes |
|---|---|---|
| Delete legacy GUI stacks (`gui/app.py`, `app_pro.py`, `main_window*`, dialogs, 10 widgets) | ✅ | Done in Phase A, -6,656 LOC. |
| Redirect `gui/__init__.launch_gui` → `workbench_app.main`; add `alprotein` console script | ✅ | Done in Phase A. |
| Split `workbench_window.py` into thin assembler + controllers | ❌ | Not done. The file is still ~2,000 lines after additions. Roadmap P0.3. |
| Spectrum widget on pyqtgraph; matplotlib only for export | 🟡 | pyqtgraph widget exists side-by-side. Matplotlib widget not yet removed. |
| Inspector, Cross-highlighting, Domains panel | 🟡 | Inspector ✅, Cross-highlighting ✅ for pigment, but Domains panel is just an inspector field, not a separate card. |
| Fit-to-experiment tab | ❌ | Phase E — not started. |

---

## Section 6 — Visible Gaps Ranked by "How much does this make it not look like the mockup?"

1. **No toolbar row.** The biggest visual difference. Adding it would convert ~10 menu/palette commands into one-click visible buttons.
2. **Mockup-style Hamiltonian heatmap.** A coloured N×N grid with click → selected pair card. Currently we have a bar plot + table.
3. **Right-sidebar cards for Live View Tools, Notebook, Domains.** Currently the right sidebar only has the Inspector.
4. **Per-exciton oscillator-strength bar plot under the spectrum.**
5. **Fit-to-Experiment tab.**
6. **Pin / Lock buttons + view-state extras on the Inspector.**
7. **"Spectra (Fast)" tab badge + retiring the matplotlib widget** so users don't see two Spectra tabs.
8. **Auto-save timer + visible "Auto-save ✓" chip + dirty marker in title bar.**
9. **History panel for Undo stack.**
10. **Help side-panel rendering full citation markdown.**

Things that are *fully done* but weren't visible because nothing was loaded in the screenshot:

- Pigment selection / cross-highlighting (needs a PDB loaded)
- Exciton scrub bar action (needs a Hamiltonian computed)
- Cancel button in status bar (only visible while running)
- Toasts (only appear on calculation errors / successes)
- Theme toggle, ⌘K (you have to invoke them)
- Tooltips with citations (have to hover)

---

## Recommendation — three slices of follow-up work

Pick whichever, all, or none.

### Slice 1 — "Make it look like the mockup" (≈1 day, low risk, mostly UI)

- Add the toolbar row (Open / Save / Recent / Undo / Redo / Run All / Stop / Export / Theme / ⌘K). 30 min.
- Add a Notebook QPlainTextEdit card to the right sidebar; round-trip through `Project.notebook`. 30 min.
- Add a Domains card with coloured swatches and "Cutoff ≥" slider in the right sidebar. 1 h.
- Add a "Live View Tools" card (color-by dropdown, μ vectors / labels toggles). 1 h.
- Dirty-marker in window title (`• modified`). 5 min.
- Make the drop zone accept `.alproj`. 10 min.
- Wire the orphaned `view_all_btn` or remove it. 5 min.
- Retire the matplotlib spectrum tab, promote the pyqtgraph one to `📊 Spectra`. 30 min.

### Slice 2 — "Build the missing high-leverage features" (≈3 days, scientific value)

- Hamiltonian heatmap with click → selected-pair card (this is a big visual + analytic upgrade).
- Per-exciton oscillator-strength bar plot under the spectrum.
- Auto-save every 60 s + dirty marker.
- Visible undo-history panel.

### Slice 3 — "Phase E + retire technical debt"

- Fit-to-experiment tab (the scipy.optimize.least_squares against loaded experimental CSV).
- Split `workbench_window.py` into controllers (the long-promised P0.3 from `ROADMAP.md`).

Tell me which slice (or which individual items) to tackle next.
