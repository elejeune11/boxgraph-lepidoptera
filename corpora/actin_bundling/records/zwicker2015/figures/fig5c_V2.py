"""Calibrated reader for Fig. 5(c) active PDE lower droplet volume trajectory V_2/w^3 vs t*k_0."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig5


def read(image, at):
    return _fig5.read_fig5ac_V(image, "V2", [6000.0, 4000.0, 2000.0])
