#!/usr/bin/env python3
"""
Minimal example: launch the Alprotein Scientific Workbench.

After ``pip install -e .`` the same thing is available as the console script
``alprotein``.
"""

from Alprotein.gui import launch_gui


if __name__ == "__main__":
    launch_gui()
