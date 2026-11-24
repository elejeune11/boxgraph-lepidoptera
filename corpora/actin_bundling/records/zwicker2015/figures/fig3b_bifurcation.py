"""Calibrated reader for Fig. 3(b) saddle-node bifurcation point (k_f^c/k_0, R_*/w) (blue dot)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig3


def read(image, at):
    return _fig3.read_fig3bd(image, at, "bifurcation")
