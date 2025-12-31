#!/usr/bin/env python3
"""
Example usage of the Alprotein GUI

This script demonstrates different ways to use the Alprotein GUI module.
"""

import sys
import os
from pathlib import Path

# Add Alprotein to path if needed
sys.path.insert(0, str(Path(__file__).parent))

def example_1_basic_launch():
    """Example 1: Basic GUI launch"""
    print("🚀 Example 1: Basic GUI Launch")
    print("-" * 30)
    
    try:
        # Import and launch GUI
        from Alprotein.gui import launch_gui
        
        print("Launching Alprotein GUI...")
        gui = launch_gui()
        
        print("✅ GUI launched successfully!")
        return gui
        
    except ImportError as e:
        print(f"❌ Import error: {e}")
        print("Make sure PyQt5 is installed: pip install PyQt5")
        return None
    except Exception as e:
        print(f"❌ Error: {e}")
        return None

def example_2_convenience_function():
    """Example 2: Using convenience function"""
    print("🎯 Example 2: Convenience Function")
    print("-" * 30)
    
    try:
        # Import and use convenience function
        from Alprotein import gui
        
        print("Launching GUI with convenience function...")
        app = gui()
        
        if app:
            print("✅ GUI launched successfully!")
        else:
            print("❌ GUI dependencies not available")
        
        return app
        
    except ImportError as e:
        print(f"❌ Import error: {e}")
        return None

def example_3_programmatic_control():
    """Example 3: Programmatic control of GUI"""
    print("🔧 Example 3: Programmatic Control")
    print("-" * 30)
    
    try:
        from PyQt5.QtWidgets import QApplication
        from Alprotein.gui import AlproteinGUI
        
        # Create application
        app = QApplication.instance()
        if app is None:
            app = QApplication(sys.argv)
        
        # Create GUI
        gui = AlproteinGUI()
        
        # Connect to signals for monitoring
        def on_structure_loaded(success, message):
            print(f"📁 Structure loaded: {success} - {message}")
        
        def on_calculation_finished(calc_type, success, message):
            print(f"⚡ Calculation finished: {calc_type} - {success} - {message}")
        
        gui.main_window.structure_loaded.connect(on_structure_loaded)
        gui.main_window.calculation_finished.connect(on_calculation_finished)
        
        # Show GUI
        gui.show()
        
        print("✅ GUI created with signal monitoring")
        print("Load a PDB file to see the signals in action")
        
        return gui
        
    except ImportError as e:
        print(f"❌ Import error: {e}")
        return None

def example_4_check_availability():
    """Example 4: Check GUI availability"""
    print("🔍 Example 4: Check GUI Availability")
    print("-" * 30)
    
    try:
        # Check if GUI is available
        import Alprotein
        
        if hasattr(Alprotein, 'launch_gui'):
            print("✅ GUI module available")
            
            # Check for web engine
            try:
                from PyQt5.QtWebEngineWidgets import QWebEngineView
                print("✅ Full 3D visualization available (QtWebEngine)")
            except ImportError:
                print("⚠️  Fallback 3D visualization (QtWebEngine not available)")
            
            # Check other dependencies
            try:
                import matplotlib
                print("✅ Plotting capabilities available")
            except ImportError:
                print("❌ Matplotlib not available")
            
            try:
                import pandas
                print("✅ Enhanced data handling available")
            except ImportError:
                print("💡 Pandas not available (optional)")
                
        else:
            print("❌ GUI module not available")
            print("Install dependencies: pip install PyQt5 matplotlib")
            
    except ImportError:
        print("❌ Alprotein package not found")

def show_system_info():
    """Show system information relevant to GUI"""
    print("💻 System Information")
    print("-" * 30)
    
    import platform
    print(f"Python: {sys.version}")
    print(f"Platform: {platform.platform()}")
    
    # Check GUI dependencies
    deps = [
        ('PyQt5', 'PyQt5'),
        ('QtWebEngine', 'PyQt5.QtWebEngineWidgets'),
        ('matplotlib', 'matplotlib'),
        ('numpy', 'numpy'),
        ('biopython', 'Bio'),
        ('pandas', 'pandas'),
        ('seaborn', 'seaborn')
    ]
    
    print("\nDependencies:")
    for name, module in deps:
        try:
            __import__(module)
            print(f"  ✅ {name}")
        except ImportError:
            print(f"  ❌ {name}")

def main():
    """Main function to run examples"""
    print("🧬 Alprotein GUI Examples")
    print("=" * 50)
    
    # Show system info
    show_system_info()
    print()
    
    # Check availability first
    example_4_check_availability()
    print()
    
    # Ask user which example to run
    print("Available examples:")
    print("1. Basic GUI launch")
    print("2. Convenience function")
    print("3. Programmatic control")
    print("0. Exit")
    
    try:
        choice = input("\nEnter your choice (0-3): ").strip()
        
        if choice == "1":
            return example_1_basic_launch()
        elif choice == "2":
            return example_2_convenience_function()
        elif choice == "3":
            return example_3_programmatic_control()
        elif choice == "0":
            print("👋 Goodbye!")
            return None
        else:
            print("❌ Invalid choice")
            return None
            
    except KeyboardInterrupt:
        print("\n👋 Goodbye!")
        return None

if __name__ == "__main__":
    main()
