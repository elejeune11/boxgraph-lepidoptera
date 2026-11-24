"""Calibrated reader for Fig. 9(a) autocatalytic PDE coarsening rates for small system (orange discs, V_s/w^3 = 1.5e5)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig9


def read(image, at):
    return _fig9.read_fig9a(image, at, "pde_small")
