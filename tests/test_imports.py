#!/usr/bin/env python3
"""Test script to debug import issues"""

import sys
import os
from pathlib import Path

print("Python version:", sys.version)
print("Current working directory:", os.getcwd())

# Add the project root to path
project_root = "/Users/mohamed/Documents/Research/Projects/2025/Alprotein-Alpha"
sys.path.insert(0, project_root)
print(f"Added {project_root} to Python path")

print("\nPython path:")
for i, path in enumerate(sys.path[:5]):  # Show first 5 paths
    print(f"  {i}: {path}")

print("\nTesting imports...")

# Test 1: Import Alprotein package
try:
    import Alprotein
    print("✅ Alprotein package imported successfully")
    print(f"   Location: {Alprotein.__file__}")
    print(f"   Package contents: {dir(Alprotein)}")
except ImportError as e:
    print(f"❌ Failed to import Alprotein: {e}")

# Test 2: Check if gui module exists
try:
    alprotein_gui_path = os.path.join(project_root, "Alprotein", "gui")
    print(f"\nChecking GUI module at: {alprotein_gui_path}")
    print(f"   Directory exists: {os.path.exists(alprotein_gui_path)}")
    print(f"   Is directory: {os.path.isdir(alprotein_gui_path)}")
    
    if os.path.exists(alprotein_gui_path):
        gui_contents = os.listdir(alprotein_gui_path)
        print(f"   Contents: {gui_contents}")
        
        # Check for __init__.py
        init_file = os.path.join(alprotein_gui_path, "__init__.py")
        print(f"   __init__.py exists: {os.path.exists(init_file)}")
except Exception as e:
    print(f"❌ Error checking GUI module: {e}")

# Test 3: Try to import gui module
try:
    from Alprotein.gui import launch_gui
    print("✅ GUI module imported successfully")
    print(f"   launch_gui function: {launch_gui}")
except ImportError as e:
    print(f"❌ Failed to import GUI module: {e}")
    import traceback
    traceback.print_exc()

print("\nTest complete!")
