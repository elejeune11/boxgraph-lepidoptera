"""Calibrated reader for Fig. 5(c) active PDE upper droplet volume trajectory V_1/w^3 vs t*k_0."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig5


def read(image, at):
    return _fig5.read_fig5ac_V(image, "V1", [6000.0, 4000.0, 2000.0])
