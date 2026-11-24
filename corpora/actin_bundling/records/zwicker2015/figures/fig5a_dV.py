"""Calibrated reader for Fig. 5(a) passive PDE volume difference Delta V/w^3 markers vs t*k_0."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig5


def read(image, at):
    return _fig5.read_fig5ac_dV(image, at, "markers")
