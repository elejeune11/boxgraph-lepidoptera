"""Calibrated reader for Fig. 3D (p. 3): overtwist Delta_theta(R*) at c_a = 2 mg/ml (red squares)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3d import read_fig3d_series


def read(image, at):
    return read_fig3d_series(image, "ca2")
