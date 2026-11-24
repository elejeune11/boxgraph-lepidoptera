"""Calibrated reader for Fig. 3(c) finite-system passive droplet growth rate dR/dt/(w*k_0) vs R/w (orange curve)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig3


def read(image, at):
    return _fig3.read_fig3ac(image, at, "passive")
