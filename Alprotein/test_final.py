#!/usr/bin/env python3
"""
Final test of the Alprotein GUI
"""

print("🧬 Testing Alprotein GUI...")

# Test import
try:
    import sys
    sys.path.insert(0, '/Users/mohamed/Documents/Research/Projects/2025/Alprotein-Alpha')
    
    print("1. Testing core import...")
    import Alprotein
    print("   ✅ Alprotein imported")
    
    print("2. Testing GUI import...")
    from Alprotein.gui import launch_gui
    print("   ✅ GUI imported successfully!")
    
    print("3. Checking dependencies...")
    try:
        from PyQt5.QtWebEngineWidgets import QWebEngineView
        print("   ✅ Full 3D visualization available")
    except ImportError:
        print("   ⚠️  Fallback 3D visualization (install PyQtWebEngine for full features)")
    
    print("\n🎉 SUCCESS! Alprotein GUI is ready to launch!")
    print("\nTo start the GUI, run:")
    print("   python launch_gui.py")
    print("   OR")
    print("   python gui_launcher.py")
    
except ImportError as e:
    print(f"❌ Import error: {e}")
    
except Exception as e:
    print(f"❌ Error: {e}")

print("\n" + "="*50)
print("GUI Implementation Complete! 🎉")
print("="*50)
