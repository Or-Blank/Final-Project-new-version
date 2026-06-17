"""
run_gui.py
==========
Double-click this file, or run:

    python run_gui.py

from inside your project folder to launch the BCR Lineage Tracer GUI.

This launcher exists because the package uses relative imports (from .loader
import ...) which only work when Python knows the folder is a package.
Running gui.py directly breaks that context — this file restores it.
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from BCR_lineage_tracer.gui import launch_gui

if __name__ == "__main__":
    launch_gui()
