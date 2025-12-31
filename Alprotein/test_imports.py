#!/usr/bin/env python3
# Quick test of GUI imports

import sys
import os
from pathlib import Path

# Add path
sys.path.insert(0, '/Users/mohamed/Documents/Research/Projects/2025/Alprotein-Alpha')

def test_imports():
    """Test all GUI imports"""
    print("🧪 Testing GUI imports...")
    
    try:
        print("1. Testing Alprotein import...")
        import Alprotein
        print("   ✅ Alprotein imported")
        
        print("2. Testing GUI module import...")
        from Alprotein.gui import launch_gui
        print("   ✅ GUI module imported")
        
        print("3. Testing main components...")
        from Alprotein.gui.app import AlproteinGUI
        from Alprotein.gui.main_window import MainWindow
        print("   ✅ Main components imported")
        
        print("4. Testing widgets...")
        from Alprotein.gui.widgets import FileLoaderWidget, CalculationPanel
        print("   ✅ Widgets imported")
        
        print("5. Testing dialogs...")
        from Alprotein.gui.dialogs import ExportDialog, SettingsDialog
        print("   ✅ Dialogs imported")
        
        print("\n🎉 All imports successful! GUI is ready to launch.")
        return True
        
    except ImportError as e:
        print(f"   ❌ Import failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    except Exception as e:
        print(f"   ❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    test_imports()
