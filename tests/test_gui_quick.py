#!/usr/bin/env python3
"""
Very simple test of GUI imports
"""

import sys
import os
from pathlib import Path

# Set up the Python path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

print("🧬 Testing Alprotein GUI imports...")

try:
    from Alprotein.gui import launch_gui
    print("✅ Successfully imported launch_gui function")
    print("🚀 Attempting to launch GUI...")
    
    # Quick test launch
    launch_gui()
    
except ImportError as e:
    print(f"❌ Import error: {e}")
    import traceback
    traceback.print_exc()
except Exception as e:
    print(f"❌ Runtime error: {e}")
    import traceback
    traceback.print_exc()
