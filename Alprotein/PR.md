```
feat: GUI Redesign and Usability Improvements

This PR introduces a significant redesign of the Alprotein GUI, focusing on improved layout, usability, and maintainability. Key changes include:

- **Interactive Layout with QSplitter**: Replaced fixed layouts with `QSplitter` in `main_window.py` to allow users to dynamically resize the sidebar and main content area, enhancing flexibility for different screen sizes and workflows.
- **Redesigned Control Sidebar**: The `ControlSidebar` has been refactored to use collapsible sections (`CollapsibleBox`) for better organization and a cleaner look. It now includes a `QScrollArea` to ensure all controls are accessible, even when content overflows the window height.
- **Restored and Enhanced Results Panel**: The "Results & Analysis" section now features a tabbed interface, explicitly restoring the "CDC Analysis" tab alongside "Site Energies." Both sections utilize a card-based layout for improved readability and visual appeal. The CDC analysis panel is now scrollable and expandable, addressing previous layout issues.
- **Improved User Feedback**: Replaced intrusive `QMessageBox` popups for successful export operations with console print statements, streamlining the user experience during repeated tasks.
- **Font Consistency and Performance**: Replaced instances of "monospace" font with "Courier" across various GUI components to resolve Qt font warnings and ensure consistent rendering.
- **Parameter Default Update**: Updated the default value for `e0a` in the calculation settings for better initial configuration.

These changes aim to provide a more modern, intuitive, and efficient user interface for Alprotein.

**Detailed Changes:**

- **`main_window.py`**:
    - Implemented `QSplitter` for main window layout, replacing `QHBoxLayout`.
    - Removed fixed width for `ControlSidebar` to enable dynamic resizing.
    - Added `init_floating_toolbar` as a placeholder for future features.
    - Set a descriptive window title.
    - Removed unused `_toggle_sidebar_collapse` method and its signal connection.
    - Replaced `QMessageBox.information` with `print` for export success messages.
- **`control_sidebar.py`**:
    - Introduced `CollapsibleBox` for organizing UI sections.
    - Integrated `QScrollArea` to manage content overflow and provide consistent spacing.
    - Corrected `PyQt.QtGui` import to `PyQt5.QtGui`.
- **`results_panel.py`**:
    - Refactored `setup_ui` to include both "Site Energies" and "CDC Analysis" tabs from initialization.
    - Implemented card-based display for results using `QFrame`.
    - Re-integrated `ContributionAnalysisWidget` and `ContributionPlotWidget` within the CDC tab.
    - Removed the redundant `add_cdc_analysis_tab` method.
    - Corrected `PyQt5.QtWidgets` import syntax.
- **`gui/dialogs/export_dialog.py`**:
    - Changed font family from "monospace" to "Courier".
- **`gui/widgets/calculation_panel.py`**:
    - Updated default value of `e0a_spinbox` to 14900.0.
    - Changed font family from "monospace" to "Courier".
- **`gui/widgets/file_loader.py`**:
    - Changed font family from "monospace" to "Courier".
- **`gui/widgets/protein_viewer.py`**:
    - Changed font family from "monospace" to "Courier".
```