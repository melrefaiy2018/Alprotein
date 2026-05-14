"""
Alprotein GUI module.

The canonical entry point is the Scientific Workbench:

    from Alprotein.gui import launch_gui
    launch_gui()

or via the installed console script ``alprotein``.
"""

from __future__ import annotations

__all__ = ["launch_gui", "main"]


def launch_gui() -> None:
    """Launch the Alprotein Scientific Workbench."""
    from .workbench_app import main as _main
    _main()


# Alias matching the entry-point convention used in pyproject.toml
main = launch_gui
