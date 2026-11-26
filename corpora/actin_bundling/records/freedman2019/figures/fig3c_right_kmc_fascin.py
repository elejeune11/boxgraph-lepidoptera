"""Calibrated reader for Fig. 3C Right — KMC fascin domain length vs polymerization rate."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _panels


def read(image, at):
    return _panels.read_fig3c_right(image, at, "fascin")
