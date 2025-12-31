#!/usr/bin/env python3
"""
Alprotein GUI Launcher

This script launches the Alprotein graphical user interface for protein
site energy calculations and visualization.

Usage:
    python gui_launcher.py [options]

Options:
    -h, --help      Show this help message
    --debug         Enable debug mode
    --test-file     Load a test PDB file on startup
"""

import sys
import os
import argparse
from pathlib import Path

# Add the Alprotein package to the Python path
current_dir = Path(__file__).parent
alprotein_dir = current_dir.parent
sys.path.insert(0, str(alprotein_dir))

def check_dependencies():
    """Check if all required dependencies are available"""
    missing_deps = []
    
    # Required packages
    required_packages = [
        ('PyQt5', 'PyQt5.QtWidgets'),
        ('matplotlib', 'matplotlib.pyplot'),
        ('numpy', 'numpy'),
        ('Bio', 'Bio.PDB'),
    ]
    
    # Optional web engine for 3D visualization
    webengine_packages = [
        ('PyQtWebEngine', 'PyQt5.QtWebEngineWidgets'),
    ]
    
    for package_name, import_name in required_packages:
        try:
            __import__(import_name)
        except ImportError:
            missing_deps.append(package_name)
    
    # Check for web engine (optional but provides better 3D visualization)
    missing_webengine = []
    for package_name, import_name in webengine_packages:
        try:
            __import__(import_name)
        except ImportError:
            missing_webengine.append(package_name)
    
    # Optional but recommended packages
    optional_packages = [
        ('pandas', 'pandas'),
        ('seaborn', 'seaborn'),
    ]
    
    missing_optional = []
    for package_name, import_name in optional_packages:
        try:
            __import__(import_name)
        except ImportError:
            missing_optional.append(package_name)
    
    if missing_deps:
        print("❌ Missing required dependencies:")
        for dep in missing_deps:
            print(f"   • {dep}")
        print("\nPlease install them using:")
        print(f"   pip install {' '.join(missing_deps)}")
        return False
    
    if missing_webengine:
        print("⚠️  Web engine not available (interactive 3D visualization disabled):")
        for dep in missing_webengine:
            print(f"   • {dep}")
        print("\nFor full 3D visualization, install:")
        print(f"   pip install {' '.join(missing_webengine)}")
        print("\nThe GUI will still work with fallback visualization.\n")
    
    if missing_optional:
        print("💡 Optional dependencies (for enhanced functionality):")
        for dep in missing_optional:
            print(f"   • {dep}")
        print(f"\nInstall them with:")
        print(f"   pip install {' '.join(missing_optional)}")
    
    return True


def main():
    """Main entry point for the GUI launcher"""
    parser = argparse.ArgumentParser(
        description="Launch Alprotein GUI for protein site energy calculations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python gui_launcher.py                    # Launch GUI normally
    python gui_launcher.py --debug            # Launch with debug output
    python gui_launcher.py --test-file        # Launch with test file
        """
    )
    
    parser.add_argument('--debug', action='store_true',
                       help='Enable debug mode with verbose output')
    parser.add_argument('--test-file', metavar='PDB_FILE',
                       help='Load specified PDB file on startup')
    parser.add_argument('--version', action='version', version='Alprotein GUI 1.0')
    
    args = parser.parse_args()
    
    # Set up logging level
    if args.debug:
        import logging
        logging.basicConfig(level=logging.DEBUG)
        print("🐛 Debug mode enabled")
    
    print("🧬 Alprotein - Protein Site Energy Calculator")
    print("=" * 50)
    
    # Check dependencies
    print("📦 Checking dependencies...")
    if not check_dependencies():
        sys.exit(1)
    
    print("✅ All required dependencies found")
    
    # Import and launch GUI
    try:
        print("🚀 Starting Alprotein GUI...")
        
        # Import the GUI module
        from Alprotein.gui import launch_gui
        
        # Launch the application
        gui = launch_gui()
        
        # Load test file if specified
        if args.test_file:
            if os.path.exists(args.test_file):
                print(f"📁 Loading test file: {args.test_file}")
                gui.main_window.load_pdb_file_path(args.test_file)
            else:
                print(f"❌ Test file not found: {args.test_file}")
        
        print("✅ GUI launched successfully")
        
    except ImportError as e:
        print(f"❌ Failed to import Alprotein GUI modules: {e}")
        print("\nPlease ensure that:")
        print("1. You're running this script from the correct directory")
        print("2. The Alprotein package is properly installed")
        print("3. All dependencies are installed")
        sys.exit(1)
        
    except Exception as e:
        print(f"❌ Failed to launch GUI: {e}")
        if args.debug:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
