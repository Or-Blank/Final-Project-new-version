import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from BCR_lineage_tracer.gui import launch_gui

if __name__ == "__main__":
    launch_gui()
