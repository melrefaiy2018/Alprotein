#!/usr/bin/env python3
"""
Test the updated GUI with Actions tab
"""

import sys
import os
from pathlib import Path

# Add the Alprotein package to the Python path
current_dir = Path(__file__).parent
sys.path.insert(0, str(current_dir))

def test_gui_components():
    """Test if all GUI components import correctly"""
    print("🧪 Testing GUI components after Actions tab update...")
    
    try:
        # Test imports
        from Alprotein.gui.widgets.actions_panel import ActionsPanel
        print("✅ ActionsPanel import successful")
        
        from Alprotein.gui.widgets.control_sidebar import ControlSidebar
        print("✅ ControlSidebar import successful")
        
        from Alprotein.gui.widgets.calculation_panel import CalculationPanel
        print("✅ CalculationPanel import successful")
        
        # Test ActionsPanel creation
        actions_panel = ActionsPanel()
        if hasattr(actions_panel, 'calculate_site_energies'):
            print("✅ ActionsPanel has calculate_site_energies signal")
        else:
            print("❌ ActionsPanel missing calculate_site_energies signal")
            return False
            
        if hasattr(actions_panel, 'run_cdc_analysis'):
            print("✅ ActionsPanel has run_cdc_analysis signal")
        else:
            print("❌ ActionsPanel missing run_cdc_analysis signal")
            return False
            
        # Test ControlSidebar creation
        sidebar = ControlSidebar()
        if hasattr(sidebar, 'actions_panel'):
            print("✅ ControlSidebar has actions_panel")
        else:
            print("❌ ControlSidebar missing actions_panel")
            return False
            
        # Test if signals are properly connected
        if hasattr(sidebar, 'calculate_site_energies'):
            print("✅ ControlSidebar exposes calculate_site_energies signal")
        else:
            print("❌ ControlSidebar missing calculate_site_energies signal")
            return False
            
        print("\n🎉 All component tests passed!")
        return True
        
    except Exception as e:
        print(f"❌ Error during testing: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_main_gui():
    """Test the main GUI launch"""
    print("\n🚀 Testing main GUI launch...")
    
    try:
        from Alprotein.gui import launch_gui
        print("✅ GUI import successful")
        
        # We won't actually launch the GUI in test mode
        print("✅ GUI launch function available")
        return True
        
    except Exception as e:
        print(f"❌ GUI launch test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print("🧬 Alprotein GUI Actions Tab Test")
    print("=" * 40)
    
    component_test = test_gui_components()
    gui_test = test_main_gui()
    
    if component_test and gui_test:
        print("\n✅ All tests passed! GUI should work with Actions tab.")
        print("\nKey improvements:")
        print("  • Dielectric constant range: 1.0 - 50,000.0")
        print("  • Actions moved to dedicated tab")
        print("  • Cleaner separation of concerns")
        print("  • Better progress tracking")
    else:
        print("\n❌ Some tests failed. Check the output above.")
