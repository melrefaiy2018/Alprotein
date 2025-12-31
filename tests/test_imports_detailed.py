#!/usr/bin/env python3
"""Quick import test for Alprotein GUI"""

import sys
import os
from pathlib import Path

# Set up the Python path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

print(f"Project root: {project_root}")
print(f"Python version: {sys.version}")

def test_import(module_name, description=""):
    """Test importing a module"""
    try:
        module = __import__(module_name, fromlist=[''])
        print(f"✅ {module_name} {description}")
        return True
    except ImportError as e:
        print(f"❌ {module_name} {description}: {e}")
        return False
    except Exception as e:
        print(f"⚠️  {module_name} {description}: {e}")
        return False

# Test basic imports
print("\n=== Basic Imports ===")
test_import("PyQt5", "(GUI framework)")
test_import("numpy", "(numerical computing)")
test_import("matplotlib", "(plotting)")
test_import("Bio", "(biopython)")

# Test Alprotein imports
print("\n=== Alprotein Core Imports ===")
if test_import("Alprotein", "(package)"):
    test_import("Alprotein.core", "(core module)")
    test_import("Alprotein.core.protein_structure", "(protein structure)")
    test_import("Alprotein.core.pigment_system", "(pigment system)")
    test_import("Alprotein.calculators", "(calculators module)")
    test_import("Alprotein.calculators.site_energy_calculator", "(site energy calculator)")

# Test GUI imports step by step
print("\n=== Alprotein GUI Imports ===")
if test_import("Alprotein.gui", "(GUI package)"):
    test_import("Alprotein.gui.widgets", "(widgets module)")
    test_import("Alprotein.gui.dialogs", "(dialogs module)")
    test_import("Alprotein.gui.app", "(main app)")
    test_import("Alprotein.gui.main_window", "(main window)")

# Test specific widgets
print("\n=== Widget Imports ===")
test_import("Alprotein.gui.widgets.file_loader", "(file loader)")
test_import("Alprotein.gui.widgets.calculation_panel", "(calculation panel)")
test_import("Alprotein.gui.widgets.results_panel", "(results panel)")
test_import("Alprotein.gui.widgets.protein_viewer", "(protein viewer)")
test_import("Alprotein.gui.widgets.control_sidebar", "(control sidebar)")
test_import("Alprotein.gui.widgets.actions_panel", "(actions panel)")
test_import("Alprotein.gui.widgets.control_panel_3d", "(3D control panel)")

# Test launch function
print("\n=== Launch Function ===")
try:
    from Alprotein.gui import launch_gui
    print("✅ launch_gui function imported successfully")
    print(f"   Function: {launch_gui}")
except ImportError as e:
    print(f"❌ launch_gui function: {e}")
    import traceback
    traceback.print_exc()

print("\n=== Test Complete ===")
