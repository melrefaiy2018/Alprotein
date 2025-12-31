#!/usr/bin/env python3

# Simple test script to run the GUI launcher
import subprocess
import sys
import os

# Get the current directory 
current_dir = os.path.dirname(os.path.abspath(__file__))
launcher_path = os.path.join(current_dir, "gui_launcher.py")

print("Testing Alprotein GUI Launcher...")
print("=" * 40)

# Run the launcher
try:
    result = subprocess.run([sys.executable, launcher_path], 
                          cwd=current_dir, 
                          capture_output=True, 
                          text=True, 
                          timeout=30)
    
    print("STDOUT:")
    print(result.stdout)
    
    if result.stderr:
        print("STDERR:")
        print(result.stderr)
    
    print(f"Return code: {result.returncode}")
    
except subprocess.TimeoutExpired:
    print("⏰ GUI launcher timed out (this is expected if GUI opened successfully)")
except Exception as e:
    print(f"❌ Error running launcher: {e}")
