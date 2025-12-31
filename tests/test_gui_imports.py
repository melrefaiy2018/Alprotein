#!/usr/bin/env python3
"""
Test script to verify all GUI imports work correctly
"""

print("Testing GUI imports...")

try:
    print("1. Importing workflow_panel...")
    from Alprotein.gui.widgets.workflow_panel import WorkflowPanel
    print("   ✓ WorkflowPanel imported successfully")

    print("2. Importing parameters_panel...")
    from Alprotein.gui.widgets.parameters_panel import ParametersPanel
    print("   ✓ ParametersPanel imported successfully")

    print("3. Importing spectrum_widget...")
    from Alprotein.gui.widgets.spectrum_widget import SpectrumPlotWidget
    print("   ✓ SpectrumPlotWidget imported successfully")

    print("4. Importing hamiltonian_widget...")
    from Alprotein.gui.widgets.hamiltonian_widget import HamiltonianWidget
    print("   ✓ HamiltonianWidget imported successfully")

    print("5. Importing enhanced_results_panel...")
    from Alprotein.gui.widgets.enhanced_results_panel import EnhancedResultsPanel
    print("   ✓ EnhancedResultsPanel imported successfully")

    print("6. Importing main_window_pro...")
    from Alprotein.gui.main_window_pro import MainWindowPro
    print("   ✓ MainWindowPro imported successfully")

    print("7. Importing app_pro...")
    from Alprotein.gui.app_pro import AlproteinGUIPro
    print("   ✓ AlproteinGUIPro imported successfully")

    print("\n✅ All imports successful!")
    print("\nThe GUI should be ready to launch with:")
    print("  python launch_pro_gui.py")

except Exception as e:
    print(f"\n❌ Import failed: {e}")
    import traceback
    traceback.print_exc()
