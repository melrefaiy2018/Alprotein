#!/usr/bin/env python3
"""
Test script to verify all Scientific Workbench imports
"""

def test_imports():
    """Test all workbench component imports"""
    print("Testing Workbench imports...\n")

    errors = []

    # Test individual widgets
    try:
        from Alprotein.gui.widgets.tools_panel import ToolsPanel
        print("✓ ToolsPanel imported successfully")
    except Exception as e:
        print(f"✗ ToolsPanel import failed: {e}")
        errors.append(f"ToolsPanel: {e}")

    try:
        from Alprotein.gui.widgets.analysis_dashboard import AnalysisDashboard
        print("✓ AnalysisDashboard imported successfully")
    except Exception as e:
        print(f"✗ AnalysisDashboard import failed: {e}")
        errors.append(f"AnalysisDashboard: {e}")

    try:
        from Alprotein.gui.widgets.data_table_widget import DataTableWidget
        print("✓ DataTableWidget imported successfully")
    except Exception as e:
        print(f"✗ DataTableWidget import failed: {e}")
        errors.append(f"DataTableWidget: {e}")

    # Test main window
    try:
        from Alprotein.gui.workbench_window import ScientificWorkbenchWindow
        print("✓ ScientificWorkbenchWindow imported successfully")
    except Exception as e:
        print(f"✗ ScientificWorkbenchWindow import failed: {e}")
        errors.append(f"ScientificWorkbenchWindow: {e}")

    # Test app
    try:
        from Alprotein.gui.workbench_app import ScientificWorkbenchApp
        print("✓ ScientificWorkbenchApp imported successfully")
    except Exception as e:
        print(f"✗ ScientificWorkbenchApp import failed: {e}")
        errors.append(f"ScientificWorkbenchApp: {e}")

    # Summary
    print("\n" + "="*60)
    if not errors:
        print("✅ All imports successful!")
        print("\nYou can now launch the workbench with:")
        print("    python launch_workbench.py")
        return True
    else:
        print(f"❌ {len(errors)} import(s) failed:")
        for error in errors:
            print(f"  - {error}")
        return False


if __name__ == "__main__":
    import sys
    success = test_imports()
    sys.exit(0 if success else 1)
