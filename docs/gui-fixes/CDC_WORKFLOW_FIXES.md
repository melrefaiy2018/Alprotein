# CDC Analysis Workflow Fixes

## Issues Addressed

### 1. **CDC Analysis Not Working**
**Problem:** CDC analysis button was enabled but the workflow was confusing
**Solution:** 
- Added prerequisite checking for site energies
- CDC button now disabled until site energies are calculated
- Clear workflow: Load → Calculate Site Energies → Run CDC Analysis

### 2. **Confusing File Dialog**
**Problem:** CDC analysis always asked for output directory, confusing users
**Solution:**
- Made file output optional
- Added choice dialog: "visualization only" vs "detailed reports"
- Users can run CDC analysis without saving files

### 3. **macOS NSOpenPanel Warning**
**Problem:** `2025-06-29 23:20:21.233 python[64288:18940811] The class 'NSOpenPanel' overrides the method identifier`
**Solution:**
- Added `QT_MAC_WANTS_LAYER=1` environment variable
- Using `QFileDialog.DontUseNativeDialog` option on macOS
- Applied to all file dialogs in the application

## Code Changes Made

### `/gui/main_window.py`
```python
# Enhanced CDC analysis workflow
def run_cdc_analysis(self):
    # Check if site energies calculated first
    if 'site_energies' not in self.results:
        # Offer to calculate site energies
        reply = QMessageBox.question(...)
    
    # Optional file output
    reply = QMessageBox.question(
        "Do you want to generate detailed file reports?"
    )
    
    # Use non-native dialog on macOS
    options = QFileDialog.DontUseNativeDialog if sys.platform == 'darwin' else QFileDialog.Options()
```

### `/gui/widgets/calculation_panel.py`
```python
# Added site energies state tracking
def __init__(self):
    self.site_energies_calculated = False

# Dynamic button state management
def update_button_states(self):
    self.cdc_analysis_button.setEnabled(
        self.system_loaded and self.site_energies_calculated
    )

# Status message updates
def set_site_energies_calculated(self, calculated: bool):
    self.site_energies_calculated = calculated
    # Update info label dynamically
```

### `/gui/app.py`
```python
# macOS warning suppression
def main():
    if sys.platform == 'darwin':
        os.environ['QT_MAC_WANTS_LAYER'] = '1'
```

### `/gui/widgets/file_loader.py`
```python
# Non-native dialogs on macOS
def browse_files(self):
    options = QFileDialog.DontUseNativeDialog if sys.platform == 'darwin' else QFileDialog.Options()
    file_path, _ = QFileDialog.getOpenFileName(..., options=options)
```

## New Workflow

### Step 1: Load Structure
- User loads PDB file
- System validates structure
- Site energies button becomes enabled
- CDC analysis button remains disabled

### Step 2: Calculate Site Energies
- User clicks "Calculate Site Energies"
- Calculation runs in background
- Upon completion:
  - CDC analysis button becomes enabled
  - Status message updates to "✅ Site energies calculated - CDC analysis ready"

### Step 3: Run CDC Analysis
- User clicks "Run CDC Analysis"
- System presents options:
  - **"Visualization only"**: Runs analysis for display in GUI
  - **"Detailed reports"**: Asks for output directory and saves files
  - **"Cancel"**: Cancels operation
- Analysis runs based on user choice

## User Experience Improvements

### Clear Prerequisites
- CDC button disabled until site energies calculated
- Helpful status messages guide user through workflow
- Option to automatically calculate site energies if missing

### Flexible Output Options
- Can run CDC analysis just for visualization
- Optional detailed file reports
- No forced directory selection

### Better Error Handling
- Clear messages when prerequisites not met
- Graceful handling of user cancellation
- Proper state management across workflow

### macOS Compatibility
- Eliminated NSOpenPanel warnings
- Consistent dialog behavior across platforms
- Better native integration

## Testing Instructions

1. **Load a PDB file**
   - Verify CDC button is disabled
   - Check status message shows site energies needed

2. **Calculate site energies**
   - Click "Calculate Site Energies"
   - Verify CDC button becomes enabled after completion
   - Check status message updates

3. **Run CDC analysis**
   - Click "Run CDC Analysis"
   - Verify option dialog appears
   - Test both "visualization only" and "detailed reports" paths

4. **macOS testing**
   - Run on macOS system
   - Verify no NSOpenPanel warnings in console
   - Check all file dialogs work properly

## Status: ✅ Complete

All CDC analysis workflow issues have been resolved:
- ✅ CDC analysis working properly
- ✅ Clear prerequisite workflow
- ✅ Optional file output
- ✅ macOS warnings eliminated
- ✅ Better user experience
