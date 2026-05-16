"""
Scientific Workbench Application

Main application for the Alprotein Scientific Workbench with:
- Professional menu bar
- Status bar with progress tracking
- Keyboard shortcuts
- Help system
"""

import sys
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QAction, QActionGroup, QMessageBox,
    QProgressBar, QLabel, QFileDialog, QMenu, QPushButton, QShortcut,
    QSizePolicy, QToolBar, QToolButton, QWidget
)
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QKeySequence, QIcon, QFont
from pathlib import Path

from Alprotein.gui.workbench_window import ScientificWorkbenchWindow
from Alprotein.gui.styles import GLOBAL_STYLESHEET, CHART_STYLE  # noqa: F401 — kept for compat
from Alprotein.gui.theme import Theme, get_theme, light_theme, dark_theme
from Alprotein.gui.widgets.toast import ToastManager
from Alprotein.gui.widgets.command_palette import CommandPalette
from Alprotein.gui.project import Project, PROJECT_SUFFIX
from Alprotein.gui.recent_files import RecentFiles
from Alprotein.gui.undo import UndoController
import matplotlib.pyplot as plt

class ScientificWorkbenchApp(QMainWindow):
    """
    Main application window for Scientific Workbench
    """

    def __init__(self):
        # IMPORTANT: setUnifiedTitleAndToolBarOnMac must be called *before*
        # the QToolBar is added; otherwise the dark macOS title bar seam
        # remains. Set it before super-init does any window construction.
        super().__init__()
        self.setUnifiedTitleAndToolBarOnMac(True)
        self.setWindowTitle("Alprotein Scientific Workbench")
        self.setGeometry(50, 50, 1600, 1000)

        # Theme state — light by default, switchable from View menu.
        self._theme: Theme = light_theme

        # Central widget
        self.workbench = ScientificWorkbenchWindow()
        self.setCentralWidget(self.workbench)

        # Non-blocking notifications (anchored bottom-right of this window).
        self.toasts = ToastManager(self, theme=self._theme)

        # ⌘K command palette
        self.palette = CommandPalette(self, theme=self._theme)

        # Project + recents + undo
        self.recent = RecentFiles()
        self.undo = UndoController(self)
        self._current_project_path: "Path | None" = None
        self._project_dirty: bool = False
        self._project_display_name: str = ""

        # Setup UI components
        self.setup_menu_bar()
        self.setup_toolbar()
        self.setup_status_bar()
        self.setup_shortcuts()
        self.connect_signals()
        self.register_palette_commands()

        # Apply styling
        self.apply_styling()

    def setup_menu_bar(self):
        """Create professional menu bar (styling driven by the active theme)."""
        menubar = self.menuBar()

        # File Menu
        file_menu = menubar.addMenu("File")

        self.open_action = QAction("Open PDB...", self)
        self.open_action.setShortcut(QKeySequence.Open)
        self.open_action.setStatusTip("Open a PDB file")
        self.open_action.triggered.connect(self.workbench.on_open_project)
        file_menu.addAction(self.open_action)

        self.open_project_action = QAction("Open Project…", self)
        self.open_project_action.setShortcut("Ctrl+Shift+O")
        self.open_project_action.setStatusTip(f"Open a saved {PROJECT_SUFFIX} project")
        self.open_project_action.triggered.connect(self.on_open_project_file)
        file_menu.addAction(self.open_project_action)

        # Open Recent submenu — populated on display.
        self.recent_menu = file_menu.addMenu("Open Recent")
        self.recent_menu.aboutToShow.connect(self._populate_recent_menu)

        file_menu.addSeparator()

        self.save_project_action = QAction("Save Project", self)
        self.save_project_action.setShortcut("Ctrl+S")
        self.save_project_action.setStatusTip("Save current session as .alproj")
        self.save_project_action.triggered.connect(self.on_save_project)
        file_menu.addAction(self.save_project_action)

        self.save_project_as_action = QAction("Save Project As…", self)
        self.save_project_as_action.setShortcut("Ctrl+Shift+S")
        self.save_project_as_action.setStatusTip("Save current session to a new file")
        self.save_project_as_action.triggered.connect(self.on_save_project_as)
        file_menu.addAction(self.save_project_as_action)

        file_menu.addSeparator()

        self.export_action = QAction("Export Data...", self)
        self.export_action.setShortcut("Ctrl+E")
        self.export_action.setStatusTip("Export calculation results")
        self.export_action.triggered.connect(self.on_export_data)
        file_menu.addAction(self.export_action)

        self.export_image_action = QAction("Export Plots...", self)
        self.export_image_action.setShortcut("Ctrl+Shift+E")
        self.export_image_action.setStatusTip("Export plots as images")
        self.export_image_action.triggered.connect(self.on_export_plots)
        file_menu.addAction(self.export_image_action)

        file_menu.addSeparator()

        self.close_action = QAction("Close Project", self)
        self.close_action.setShortcut("Ctrl+W")
        self.close_action.setStatusTip("Close current project")
        self.close_action.triggered.connect(self.on_close_project)
        file_menu.addAction(self.close_action)

        file_menu.addSeparator()

        self.exit_action = QAction("Exit", self)
        self.exit_action.setShortcut(QKeySequence.Quit)
        self.exit_action.setStatusTip("Exit application")
        self.exit_action.triggered.connect(self.close)
        file_menu.addAction(self.exit_action)

        # Edit Menu
        edit_menu = menubar.addMenu("Edit")

        self.undo_action = self.undo.make_undo_action(self)
        edit_menu.addAction(self.undo_action)

        self.redo_action = self.undo.make_redo_action(self)
        edit_menu.addAction(self.redo_action)

        edit_menu.addSeparator()

        self.settings_action = QAction("Settings...", self)
        self.settings_action.setShortcut("Ctrl+,")
        self.settings_action.setStatusTip("Open settings")
        self.settings_action.triggered.connect(self.on_settings)
        edit_menu.addAction(self.settings_action)

        # Structure Menu
        structure_menu = menubar.addMenu("Structure")

        self.properties_action = QAction("Properties...", self)
        self.properties_action.setShortcut("Ctrl+I")
        self.properties_action.setStatusTip("View structure properties")
        self.properties_action.triggered.connect(self.workbench.on_show_properties)
        structure_menu.addAction(self.properties_action)

        structure_menu.addSeparator()

        self.validate_action = QAction("Validate Structure", self)
        self.validate_action.setShortcut("Ctrl+Shift+V")
        self.validate_action.setStatusTip("Validate current structure")
        self.validate_action.triggered.connect(self.on_validate_structure)
        structure_menu.addAction(self.validate_action)

        # Calculate Menu
        calc_menu = menubar.addMenu("Calculate")

        self.calc_site_energy_action = QAction("Site Energies", self)
        self.calc_site_energy_action.setShortcut("Ctrl+1")
        self.calc_site_energy_action.setStatusTip("Calculate site energies")
        self.calc_site_energy_action.triggered.connect(
            lambda: self.workbench.on_run_calculation("site_energies")
        )
        calc_menu.addAction(self.calc_site_energy_action)

        self.calc_hamiltonian_action = QAction("Hamiltonian", self)
        self.calc_hamiltonian_action.setShortcut("Ctrl+2")
        self.calc_hamiltonian_action.setStatusTip("Construct Hamiltonian matrix")
        self.calc_hamiltonian_action.triggered.connect(
            lambda: self.workbench.on_run_calculation("hamiltonian")
        )
        calc_menu.addAction(self.calc_hamiltonian_action)

        self.calc_spectrum_action = QAction("Absorption Spectrum", self)
        self.calc_spectrum_action.setShortcut("Ctrl+3")
        self.calc_spectrum_action.setStatusTip("Calculate absorption spectrum")
        self.calc_spectrum_action.triggered.connect(
            lambda: self.workbench.on_run_calculation("spectrum")
        )
        calc_menu.addAction(self.calc_spectrum_action)

        calc_menu.addSeparator()

        self.calc_all_action = QAction("Run All Calculations", self)
        self.calc_all_action.setShortcut("Ctrl+Shift+R")
        self.calc_all_action.setStatusTip("Run complete workflow")
        self.calc_all_action.triggered.connect(self.on_run_all_calculations)
        calc_menu.addAction(self.calc_all_action)

        # Analyze Menu
        analyze_menu = menubar.addMenu("Analyze")

        self.view_site_energies_action = QAction("Site Energy Analysis", self)
        self.view_site_energies_action.setShortcut("Ctrl+Shift+1")
        self.view_site_energies_action.setStatusTip("View detailed site energy analysis")
        self.view_site_energies_action.triggered.connect(
            lambda: self.workbench.on_view_details("site_energies")
        )
        analyze_menu.addAction(self.view_site_energies_action)

        self.view_hamiltonian_action = QAction("Hamiltonian Analysis", self)
        self.view_hamiltonian_action.setShortcut("Ctrl+Shift+2")
        self.view_hamiltonian_action.setStatusTip("View detailed Hamiltonian analysis")
        self.view_hamiltonian_action.triggered.connect(
            lambda: self.workbench.on_view_details("hamiltonian")
        )
        analyze_menu.addAction(self.view_hamiltonian_action)

        self.view_spectrum_action = QAction("Spectrum Analysis", self)
        self.view_spectrum_action.setShortcut("Ctrl+Shift+3")
        self.view_spectrum_action.setStatusTip("View detailed spectrum analysis")
        self.view_spectrum_action.triggered.connect(
            lambda: self.workbench.on_view_details("spectrum")
        )
        analyze_menu.addAction(self.view_spectrum_action)

        # Visualize Menu
        viz_menu = menubar.addMenu("Visualize")

        self.reset_view_action = QAction("Reset 3D View", self)
        self.reset_view_action.setShortcut("Ctrl+R")
        self.reset_view_action.setStatusTip("Reset 3D viewer to default view")
        self.reset_view_action.triggered.connect(self.on_reset_view)
        viz_menu.addAction(self.reset_view_action)

        viz_menu.addSeparator()

        # Tab switching actions
        for shortcut, label, fragment in (
            ("Ctrl+Shift+1", "Go to 3D Structure", "3D Structure"),
            ("Ctrl+Shift+2", "Go to Hamiltonian", "Hamiltonian"),
            ("Ctrl+Shift+3", "Go to Spectra", "Spectra"),
            ("Ctrl+Shift+4", "Go to Data Analysis", "Data Analysis"),
        ):
            action = QAction(label, self)
            action.setShortcut(shortcut)
            action.setStatusTip(label)
            action.triggered.connect(
                lambda _checked=False, name=fragment: self._goto_tab_by_name(name)
            )
            viz_menu.addAction(action)

        # View Menu — theme toggle, command palette
        view_menu = menubar.addMenu("View")

        theme_menu = view_menu.addMenu("Theme")
        self.theme_action_group = QActionGroup(self)
        self.theme_action_group.setExclusive(True)

        self.theme_light_action = QAction("Light", self, checkable=True)
        self.theme_light_action.setChecked(self._theme.name == "light")
        self.theme_light_action.triggered.connect(lambda: self.set_theme("light"))
        self.theme_action_group.addAction(self.theme_light_action)
        theme_menu.addAction(self.theme_light_action)

        self.theme_dark_action = QAction("Dark", self, checkable=True)
        self.theme_dark_action.setChecked(self._theme.name == "dark")
        self.theme_dark_action.triggered.connect(lambda: self.set_theme("dark"))
        self.theme_action_group.addAction(self.theme_dark_action)
        theme_menu.addAction(self.theme_dark_action)

        self.toggle_theme_action = QAction("Toggle Light/Dark", self)
        self.toggle_theme_action.setShortcut("Ctrl+Shift+T")
        self.toggle_theme_action.triggered.connect(self.toggle_theme)
        view_menu.addAction(self.toggle_theme_action)

        view_menu.addSeparator()

        self.palette_action = QAction("Command Palette…", self)
        self.palette_action.setShortcut("Ctrl+K")
        self.palette_action.setStatusTip("Search and run any action")
        self.palette_action.triggered.connect(self.palette.toggle)
        view_menu.addAction(self.palette_action)

        # Help Menu
        help_menu = menubar.addMenu("Help")

        self.quick_start_action = QAction("Quick Start Guide", self)
        self.quick_start_action.setShortcut(QKeySequence.HelpContents)
        self.quick_start_action.setStatusTip("View quick start guide")
        self.quick_start_action.triggered.connect(self.on_quick_start)
        help_menu.addAction(self.quick_start_action)

        help_menu.addSeparator()

        self.about_action = QAction("About Alprotein", self)
        self.about_action.setStatusTip("About this application")
        self.about_action.triggered.connect(self.on_about)
        help_menu.addAction(self.about_action)

    def setup_toolbar(self):
        """Build the horizontal action toolbar shown under the menu bar.

        Items here are duplicates of menu actions and palette commands —
        the toolbar is purely a visibility / discoverability win. We
        deliberately keep the icon set as Unicode glyphs so we don't pull
        in icon-theme dependencies just for the toolbar.
        """
        bar = QToolBar("Main", self)
        bar.setObjectName("workbench_toolbar")
        bar.setMovable(False)
        bar.setFloatable(False)
        bar.setIconSize(self.iconSize())
        bar.setToolButtonStyle(Qt.ToolButtonTextOnly)
        bar.setContextMenuPolicy(Qt.PreventContextMenu)
        bar.setObjectName("workbench_toolbar")  # for QSS targeting via the theme
        self._apply_toolbar_style(bar)

        # --- File group -------------------------------------------------
        open_btn = self._toolbar_action(bar, "Open", self.workbench.on_open_project,
                                         tooltip="Open a PDB file (⌘O)", primary=True)
        save_btn = self._toolbar_action(bar, "Save", self.on_save_project,
                                         tooltip="Save the current session as .alproj (⌘S)")

        # Recent ▾ — a tool button with a popup menu that we rebuild on display.
        self._toolbar_recent_btn = QToolButton(bar)
        self._toolbar_recent_btn.setText("Recent ▾")
        self._toolbar_recent_btn.setPopupMode(QToolButton.InstantPopup)
        self._toolbar_recent_btn.setToolTip("Recently opened files")
        self._toolbar_recent_menu = QMenu(self._toolbar_recent_btn)
        self._toolbar_recent_menu.aboutToShow.connect(
            lambda: self._populate_recent_into(self._toolbar_recent_menu)
        )
        self._toolbar_recent_btn.setMenu(self._toolbar_recent_menu)
        bar.addWidget(self._toolbar_recent_btn)

        bar.addSeparator()

        # --- Edit -------------------------------------------------------
        bar.addAction(self.undo_action)
        bar.addAction(self.redo_action)

        bar.addSeparator()

        # --- Calculate --------------------------------------------------
        self._toolbar_action(bar, "Run All",
                              self.on_run_all_calculations,
                              tooltip="Run site energies → Hamiltonian → spectrum (⌘⇧R)",
                              primary=True)
        self._stop_btn = self._toolbar_action(bar, "Stop",
                                              self.cancel_current_calculation,
                                              tooltip="Cancel running calculation (⌘.)")

        bar.addSeparator()

        # --- Export / theme / palette -----------------------------------
        export_btn = QToolButton(bar)
        export_btn.setText("Export ▾")
        export_btn.setPopupMode(QToolButton.InstantPopup)
        export_menu = QMenu(export_btn)
        export_menu.addAction("Data…", self.on_export_data)
        export_menu.addAction("Plots…", self.on_export_plots)
        export_btn.setMenu(export_menu)
        bar.addWidget(export_btn)

        theme_btn = QToolButton(bar)
        theme_btn.setText("Theme ▾")
        theme_btn.setPopupMode(QToolButton.InstantPopup)
        theme_menu = QMenu(theme_btn)
        theme_menu.addAction("Light", lambda: self.set_theme("light"))
        theme_menu.addAction("Dark", lambda: self.set_theme("dark"))
        theme_btn.setMenu(theme_menu)
        bar.addWidget(theme_btn)

        # Spacer pushes the ⌘K hint to the right.
        spacer = QWidget(bar)
        spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        bar.addWidget(spacer)

        self._toolbar_action(bar, "⌘K Command palette",
                              self.palette.toggle,
                              tooltip="Search any action (⌘K)")

        self.addToolBar(Qt.TopToolBarArea, bar)
        self._toolbar = bar

    def _toolbar_action(self, bar: QToolBar, label: str, slot, *,
                         tooltip: str = "", primary: bool = False) -> QToolButton:
        """Add a plain push-button-style tool button to ``bar``."""
        btn = QToolButton(bar)
        btn.setText(label)
        if tooltip:
            btn.setToolTip(tooltip)
        btn.clicked.connect(slot)
        if primary:
            btn.setProperty("primary", "true")
        bar.addWidget(btn)
        return btn

    def _populate_recent_into(self, menu: QMenu) -> None:
        """Shared helper: fill an arbitrary QMenu with recent entries."""
        menu.clear()
        entries = self.recent.entries()
        if not entries:
            empty = menu.addAction("(no recent files)")
            empty.setEnabled(False)
            return
        for entry in entries:
            label = Path(entry.path).name
            action = menu.addAction(label)
            action.setStatusTip(entry.path)
            entry_path, entry_kind = entry.path, entry.kind
            action.triggered.connect(
                lambda _checked=False, p=entry_path, k=entry_kind: self._open_recent(p, k)
            )
        menu.addSeparator()
        clear = menu.addAction("Clear Recent")
        clear.triggered.connect(self._clear_recent)

    def setup_status_bar(self):
        """Create status bar with progress tracking and cancel button."""
        statusbar = self.statusBar()

        self.status_label = QLabel("Ready")
        self.status_label.setStyleSheet("padding: 4px 8px;")
        statusbar.addWidget(self.status_label, 1)

        self.progress_bar = QProgressBar()
        self.progress_bar.setMaximumWidth(200)
        self.progress_bar.setMaximumHeight(8)
        self.progress_bar.setVisible(False)
        statusbar.addPermanentWidget(self.progress_bar)

        self.cancel_button = QPushButton("■ Cancel")
        self.cancel_button.setProperty("class", "danger action-small")
        self.cancel_button.setVisible(False)
        self.cancel_button.setCursor(Qt.PointingHandCursor)
        self.cancel_button.clicked.connect(self.cancel_current_calculation)
        statusbar.addPermanentWidget(self.cancel_button)

        self.info_label = QLabel("")
        self.info_label.setProperty("class", "label")
        statusbar.addPermanentWidget(self.info_label)

    def setup_shortcuts(self):
        """Application-level shortcuts that don't have menu entries."""
        # Cancel running calculation
        QShortcut(QKeySequence("Ctrl+."), self, activated=self.cancel_current_calculation)

    def connect_signals(self):
        """Connect workbench signals to app."""
        self.workbench.status_message.connect(self.update_status)
        self.workbench.progress_update.connect(self.update_progress)
        self.workbench.calculation_error.connect(self.on_calculation_error)
        self.workbench.calculation_cancelled.connect(self.on_calculation_cancelled)
        self.workbench.calculation_completed.connect(self.on_calculation_completed)
        self.workbench.pdb_loaded.connect(self.on_pdb_loaded)
        # .alproj drops bubble up to the app so we can run the full project
        # restoration path instead of treating them as PDBs.
        self.workbench.project_open_requested.connect(
            lambda path: self.load_project_path(Path(path))
        )
        # Clear the undo stack when a new PDB / project is opened to avoid
        # cross-session commands that no longer make sense.
        self.workbench.pdb_loaded.connect(lambda _p: self.undo.clear())
        # Record parameter edits onto the undo stack.
        self.workbench.parameter_committed.connect(self.on_parameter_committed)
        # Mark dirty on notebook edits too.
        self.workbench.notebook_edited.connect(lambda: self._mark_dirty(True))

    def apply_styling(self):
        """Apply the current theme stylesheet plus matplotlib rcParams."""
        QApplication.instance().setStyleSheet(self._theme.qss())
        for key, value in self._theme.chart_style().items():
            plt.rcParams[key] = value
        # Re-stamp the bespoke chrome each time the theme changes — Qt's
        # widget-level stylesheets win over the application stylesheet, and
        # those bits (toolbar, status bar) need explicit token colours so
        # they don't pick up the platform's native chrome.
        if hasattr(self, "_toolbar"):
            self._apply_toolbar_style(self._toolbar)
        if self.statusBar() is not None:
            self._apply_statusbar_style()

    def _apply_toolbar_style(self, bar: QToolBar) -> None:
        """Restyle the toolbar with the current theme tokens."""
        t = self._theme
        bar.setStyleSheet(
            f"""
            QToolBar#workbench_toolbar {{
                background-color: {t.bg_panel};
                border: none;
                border-bottom: 1px solid {t.border};
                spacing: 4px;
                padding: 8px 14px;
            }}
            QToolBar#workbench_toolbar QToolButton {{
                color: {t.text_primary};
                background: transparent;
                border: 1px solid transparent;
                padding: 5px 11px;
                border-radius: 5px;
                font-size: 12px;
                font-weight: 600;
            }}
            QToolBar#workbench_toolbar QToolButton:hover {{
                background: {t.bg_subtle};
                border-color: {t.border};
            }}
            QToolBar#workbench_toolbar QToolButton[primary="true"] {{
                background: {t.accent};
                color: {t.text_inverse};
                border-color: {t.accent};
            }}
            QToolBar#workbench_toolbar QToolButton[primary="true"]:hover {{
                background: {t.accent_hover};
                border-color: {t.accent_hover};
            }}
            QToolBar#workbench_toolbar QToolButton::menu-indicator {{
                image: none;
            }}
            QToolBar#workbench_toolbar::separator {{
                background: {t.border};
                width: 1px;
                margin: 6px 4px;
            }}
            """
        )

    def _apply_statusbar_style(self) -> None:
        """Match the status bar to the theme so it doesn't read as native dark."""
        t = self._theme
        sb = self.statusBar()
        sb.setStyleSheet(
            f"""
            QStatusBar {{
                background-color: {t.bg_panel};
                color: {t.text_secondary};
                border-top: 1px solid {t.border};
            }}
            QStatusBar QLabel {{ color: {t.text_secondary}; }}
            """
        )

    # ------------------------------------------------------------------
    # Status & progress
    # ------------------------------------------------------------------

    def update_status(self, message: str):
        """Update status bar message; auto-clear after 5s for non-errors."""
        self.status_label.setText(message)
        if not message.lower().startswith("error"):
            QTimer.singleShot(5000, lambda: self.status_label.setText("Ready"))

    def update_progress(self, message: str, percentage: int):
        """Update progress bar; show the cancel button while running."""
        active = 0 < percentage < 100
        self.progress_bar.setVisible(active)
        self.cancel_button.setVisible(active)
        if active:
            self.progress_bar.setValue(percentage)
            self.status_label.setText(message)

    # ------------------------------------------------------------------
    # Calculation outcome routing
    # ------------------------------------------------------------------

    def cancel_current_calculation(self):
        """Forwarded to the workbench; toast the result."""
        if self.workbench.cancel_current_calculation():
            self.toasts.warning("Cancelling…", "Calculation will stop at the next safe point.")
        else:
            self.toasts.info("Nothing to cancel", "No calculation is running.")

    def on_calculation_error(self, calc_type: str, error_msg: str):
        """Surface worker errors as a non-blocking toast (with copyable details)."""
        first_line = error_msg.splitlines()[0] if error_msg else "Unknown error"
        self.progress_bar.setVisible(False)
        self.cancel_button.setVisible(False)
        self.toasts.error(
            f"{calc_type} failed",
            first_line[:200],
            details=error_msg,
        )

    def on_calculation_cancelled(self, calc_type: str):
        self.progress_bar.setVisible(False)
        self.cancel_button.setVisible(False)
        self.toasts.warning(f"{calc_type} cancelled", "")

    def on_calculation_completed(self, calc_type: str):
        # Hide progress when a top-level workflow step finishes.
        self.progress_bar.setVisible(False)
        self.cancel_button.setVisible(False)
        if calc_type == "workflow":
            self.toasts.success(
                "Workflow complete",
                "Site energies, Hamiltonian, spectra, and exciton distributions are ready.",
            )

    # ------------------------------------------------------------------
    # Theme
    # ------------------------------------------------------------------

    def set_theme(self, name: str):
        """Switch the active theme by name ("light" or "dark")."""
        self._theme = get_theme(name)
        self.toasts.set_theme(self._theme)
        self.apply_styling()
        if hasattr(self, "theme_light_action"):
            self.theme_light_action.setChecked(self._theme.name == "light")
            self.theme_dark_action.setChecked(self._theme.name == "dark")

    def toggle_theme(self):
        self.set_theme("dark" if self._theme.name == "light" else "light")

    # ------------------------------------------------------------------
    # Command palette
    # ------------------------------------------------------------------

    def register_palette_commands(self):
        """Register the main actions in the ⌘K command palette."""
        p, w = self.palette, self.workbench

        p.add_command("file.open", "Open PDB…", w.on_open_project,
                      group="File", shortcut="⌘O", keywords=["load", "structure"])
        p.add_command("file.open_project", "Open Project…", self.on_open_project_file,
                      group="File", shortcut="⌘⇧O", keywords=["alproj"])
        p.add_command("file.save_project", "Save Project", self.on_save_project,
                      group="File", shortcut="⌘S")
        p.add_command("file.save_project_as", "Save Project As…", self.on_save_project_as,
                      group="File", shortcut="⌘⇧S")
        p.add_command("file.export", "Export data…", self.on_export_data,
                      group="File", shortcut="⌘E")
        p.add_command("file.export_plots", "Export plots…", self.on_export_plots,
                      group="File", shortcut="⌘⇧E")
        p.add_command("file.close", "Close project", self.on_close_project,
                      group="File")

        p.add_command("edit.undo", "Undo", self.undo.stack.undo,
                      group="Edit", shortcut="⌘Z")
        p.add_command("edit.redo", "Redo", self.undo.stack.redo,
                      group="Edit", shortcut="⌘⇧Z")

        p.add_command("calc.site", "Calculate site energies",
                      lambda: w.on_run_calculation("site_energies"),
                      group="Calculate", shortcut="⌘1")
        p.add_command("calc.hamiltonian", "Build Hamiltonian",
                      lambda: w.on_run_calculation("hamiltonian"),
                      group="Calculate", shortcut="⌘2")
        p.add_command("calc.spectrum", "Calculate absorption spectrum",
                      lambda: w.on_run_calculation("spectrum"),
                      group="Calculate", shortcut="⌘3")
        p.add_command("calc.all", "Run complete workflow",
                      self.on_run_all_calculations,
                      group="Calculate", shortcut="⌘⇧R")
        p.add_command("calc.cancel", "Cancel running calculation",
                      self.cancel_current_calculation,
                      group="Calculate", shortcut="⌘.")

        p.add_command("view.theme.light", "Theme: Light",
                      lambda: self.set_theme("light"), group="View")
        p.add_command("view.theme.dark", "Theme: Dark",
                      lambda: self.set_theme("dark"), group="View")
        p.add_command("view.theme.toggle", "Toggle theme",
                      self.toggle_theme, group="View", shortcut="⌘⇧T")
        p.add_command("view.reset3d", "Reset 3D view",
                      self.on_reset_view, group="View", shortcut="⌘R")
        p.add_command("view.tab.structure", "Go to 3D Structure",
                      lambda: self._goto_tab_by_name("3D Structure"),
                      group="View", shortcut="⌘⇧1")
        p.add_command("view.tab.hamiltonian", "Go to Hamiltonian",
                      lambda: self._goto_tab_by_name("Hamiltonian"),
                      group="View", shortcut="⌘⇧2")
        p.add_command("view.tab.spectra", "Go to Spectra",
                      lambda: self._goto_tab_by_name("Spectra"),
                      group="View", shortcut="⌘⇧3")
        p.add_command("view.tab.analysis", "Go to Data & Analysis",
                      lambda: self._goto_tab_by_name("Data Analysis"),
                      group="View", shortcut="⌘⇧4")

        p.add_command("help.quickstart", "Quick start guide",
                      self.on_quick_start, group="Help", shortcut="F1")
        p.add_command("help.about", "About Alprotein", self.on_about, group="Help")

    # ------------------------------------------------------------------
    # Undo recording
    # ------------------------------------------------------------------

    def on_parameter_committed(self, key: str, old, new) -> None:
        """Record a parameter edit as an undoable command."""
        from Alprotein.gui.undo import ParameterEditCommand

        def apply(k: str, value) -> None:
            self.workbench.tools_panel.set_value_silently(k, value)

        self.undo.push(ParameterEditCommand(apply, key, old, new))
        self._mark_dirty(True)

    # ------------------------------------------------------------------
    # PDB load hook
    # ------------------------------------------------------------------

    def _goto_tab_by_name(self, fragment: str) -> None:
        """Switch to the workspace tab whose text contains ``fragment``."""
        tabs = self.workbench.workspace_tabs
        for i in range(tabs.count()):
            if fragment in tabs.tabText(i):
                tabs.setCurrentIndex(i)
                return

    def on_pdb_loaded(self, file_path: str) -> None:
        """Record the loaded PDB in recents and refresh the title."""
        path = Path(file_path)
        if path.exists():
            self.recent.add(path, kind="pdb")
        # Loading a fresh PDB invalidates the active project file.
        self._current_project_path = None
        self._project_dirty = False
        self._project_display_name = path.name
        self._refresh_title()

    def _refresh_title(self) -> None:
        """Re-render the window title from the current project state."""
        name = self._project_display_name or "Untitled"
        dirty = " • modified" if self._project_dirty else ""
        self.setWindowTitle(f"Alprotein — {name}{dirty}")

    def _mark_dirty(self, dirty: bool = True) -> None:
        if self._project_dirty == dirty:
            return
        self._project_dirty = dirty
        self._refresh_title()

    # ------------------------------------------------------------------
    # Project save / load / recents
    # ------------------------------------------------------------------

    def on_save_project(self) -> None:
        """Save to the current project path; falls back to Save As."""
        if self._current_project_path is None:
            self.on_save_project_as()
            return
        self._save_to(self._current_project_path)

    def on_save_project_as(self) -> None:
        suggested_name = "untitled"
        if self.workbench.current_file is not None:
            suggested_name = self.workbench.current_file.stem
        suggested = str(Path.home() / f"{suggested_name}{PROJECT_SUFFIX}")
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Project As",
            suggested,
            f"Alprotein Project (*{PROJECT_SUFFIX})",
        )
        if not path:
            return
        self._save_to(Path(path))

    def _save_to(self, path: Path) -> None:
        try:
            project = self.workbench.to_project()
            written = project.save(path)
        except Exception as e:  # noqa: BLE001 — surface to user
            self.toasts.error("Save failed", str(e))
            return
        self._current_project_path = written
        self._project_dirty = False
        self._project_display_name = written.name
        self.recent.add(written, kind="project")
        self._refresh_title()
        self.toasts.success("Project saved", str(written))

    def on_open_project_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Open Project",
            str(Path.home()),
            f"Alprotein Project (*{PROJECT_SUFFIX});;All files (*)",
        )
        if path:
            self.load_project_path(Path(path))

    def load_project_path(self, path: Path) -> None:
        """Open the project at ``path``."""
        try:
            project = Project.load(path)
            self.workbench.apply_project(project)
        except Exception as e:  # noqa: BLE001 — friendly toast instead of crash
            self.toasts.error("Could not open project", str(e), details=repr(e))
            return
        self._current_project_path = path
        self._project_dirty = False
        self._project_display_name = path.name
        self.recent.add(path, kind="project")
        self._refresh_title()
        self.toasts.success("Project opened", project.name or path.name)

    def _populate_recent_menu(self) -> None:
        """Rebuild the File → Open Recent submenu when it's about to display."""
        self._populate_recent_into(self.recent_menu)

    def _open_recent(self, path: str, kind: str) -> None:
        target = Path(path)
        if not target.exists():
            self.toasts.warning("File missing", f"{target.name} could not be found.")
            return
        if kind == "project":
            self.load_project_path(target)
        else:
            self.workbench.load_pdb_file(str(target))
            self.recent.add(target, kind="pdb")

    def _clear_recent(self) -> None:
        self.recent.clear()
        self.toasts.info("Recent files cleared", "")

    def on_export_data(self):
        """Export calculation data"""
        self.workbench.on_export_table()

    def on_export_plots(self):
        """Export plots as images"""
        QMessageBox.information(
            self,
            "Export Plots",
            "Click the 'Export' button on each plot in the Analysis Dashboard to export individual plots."
        )

    def on_close_project(self):
        """Close current project"""
        reply = QMessageBox.question(
            self,
            "Close Project",
            "Are you sure you want to close the current project?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )

        if reply == QMessageBox.Yes:
            self.workbench.clear_all_data()
            self.update_status("Project closed")

    def on_settings(self):
        """Open settings dialog"""
        QMessageBox.information(
            self,
            "Settings",
            "Settings are available in the Tools Panel on the right.\n\n"
            "Adjust parameters in the SETTINGS section:\n"
            "- Dielectric constants\n"
            "- Oscillator strength\n"
            "- Temperature\n"
            "- Coupling method"
        )

    def on_validate_structure(self):
        """Validate current structure"""
        if not self.workbench.pigment_system:
            QMessageBox.warning(
                self,
                "No Structure",
                "Please load a PDB file first"
            )
            return

        # Structure is validated on load, so just show info
        self.workbench.on_show_properties()

    def on_run_all_calculations(self):
        """Run complete calculation workflow"""
        if not self.workbench.pigment_system:
            QMessageBox.warning(
                self,
                "No Structure",
                "Please load a PDB file first"
            )
            return

        reply = QMessageBox.question(
            self,
            "Run All Calculations",
            "This will run the complete workflow:\n"
            "1. Site Energies\n"
            "2. Hamiltonian Construction\n"
            "3. Spectrum Calculation\n\n"
            "This may take several minutes. Continue?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )

        if reply == QMessageBox.Yes:
            # Start with site energies
            # The workflow will continue automatically via signals
            self.workbench.on_run_calculation("site_energies")

    def on_reset_view(self):
        """Reset 3D view"""
        self.workbench.protein_viewer.reset_view()
        self.update_status("3D view reset")

    def on_quick_start(self):
        """Show quick start guide"""
        guide = """
        <h2>Alprotein Scientific Workbench - Quick Start</h2>

        <h3>1. Load Structure</h3>
        <p>• File → Open PDB... (or Ctrl+O)<br>
        • Select your PDB file<br>
        • View structure in 3D viewer (first tab)</p>

        <h3>2. Configure Parameters</h3>
        <p>• Adjust settings in Tools Panel (right sidebar)<br>
        • Set dielectric constants<br>
        • Choose coupling method (TrEsp or Dipole)</p>

        <h3>3. Run Calculations</h3>
        <p>• <b>Site Energies</b>: Calculate → Site Energies (Ctrl+1)<br>
        • <b>Hamiltonian</b>: Calculate → Hamiltonian (Ctrl+2)<br>
        • <b>Spectrum</b>: Calculate → Absorption Spectrum (Ctrl+3)</p>

        <h3>4. View Results (Tabbed Interface)</h3>
        <p>• <b>Tab 1 - 3D Structure</b>: Protein visualization (Ctrl+Shift+1)<br>
        • <b>Tab 2 - Analysis Dashboard</b>: Three mini-plots (Ctrl+Shift+2)<br>
        • <b>Tab 3 - Data Table</b>: Detailed data (Ctrl+Shift+3)<br>
        • Click "Details" or "Export" in dashboard for full views</p>

        <h3>Keyboard Shortcuts</h3>
        <p><b>File Operations:</b><br>
        • <b>Ctrl+O</b>: Open PDB<br>
        • <b>Ctrl+E</b>: Export data<br>
        <b>Calculations:</b><br>
        • <b>Ctrl+1/2/3</b>: Run site energies/Hamiltonian/spectrum<br>
        <b>View Tabs:</b><br>
        • <b>Ctrl+Shift+1/2/3</b>: Switch to Structure/Dashboard/Table<br>
        <b>Other:</b><br>
        • <b>Ctrl+R</b>: Reset 3D view<br>
        • <b>F1</b>: This help</p>

        <h3>Tips</h3>
        <p>• Switch between tabs to view different aspects of your analysis<br>
        • Watch progress in status bar during calculations<br>
        • Select pigments in table to highlight in 3D viewer<br>
        • Use Ctrl+Shift+1/2/3 for quick tab switching</p>
        """

        msg = QMessageBox(self)
        msg.setWindowTitle("Quick Start Guide")
        msg.setTextFormat(Qt.RichText)
        msg.setText(guide)
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()

    def on_about(self):
        """Show about dialog"""
        about_text = """
        <h2>Alprotein Scientific Workbench</h2>
        <p><b>Version 2.0</b></p>

        <p>Professional software for pigment-protein complex spectroscopy analysis.</p>

        <h3>Features</h3>
        <ul>
        <li>Site energy calculation using CDC method</li>
        <li>Hamiltonian construction (TrEsp/Dipole coupling)</li>
        <li>Linear absorption spectrum calculation</li>
        <li>3D molecular visualization</li>
        <li>Professional data analysis and export</li>
        </ul>

        <h3>Components</h3>
        <p>• 3D Structure Viewer<br>
        • Analysis Dashboard<br>
        • Data Table<br>
        • Tools Panel</p>

        <p><small>Built with PyQt5 and Matplotlib</small></p>
        """

        msg = QMessageBox(self)
        msg.setWindowTitle("About Alprotein")
        msg.setTextFormat(Qt.RichText)
        msg.setText(about_text)
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()


def main():
    """Main application entry point"""
    import os
    # Must be set BEFORE QApplication is constructed:
    #  - AA_ShareOpenGLContexts silences the "initialized from a plugin"
    #    warning and lets QtWebEngine compose with QOpenGLWidget.
    #  - QT_MAC_WANTS_LAYER=1 puts every Qt window on a CALayer on macOS,
    #    which fixes title-bar tinting and HiDPI redraw glitches.
    from PyQt5.QtCore import QCoreApplication
    QCoreApplication.setAttribute(Qt.AA_ShareOpenGLContexts, True)
    if sys.platform == "darwin":
        os.environ.setdefault("QT_MAC_WANTS_LAYER", "1")

    app = QApplication(sys.argv)

    # Set application metadata
    app.setApplicationName("Alprotein Scientific Workbench")
    app.setOrganizationName("Alprotein")
    app.setApplicationVersion("2.0")

    # Default font — pick the OS-native UI face so Qt doesn't waste ~150 ms
    # resolving a font that isn't installed on this platform.
    if sys.platform == "darwin":
        font = QFont(".AppleSystemUIFont", 13)
    elif sys.platform.startswith("win"):
        font = QFont("Segoe UI", 10)
    else:
        font = QFont("Noto Sans", 10)
    app.setFont(font)

    # Create and show main window — apply_styling() inside the window applies
    # the active theme to QApplication and matplotlib.
    window = ScientificWorkbenchApp()
    window.show()

    # Start event loop
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
