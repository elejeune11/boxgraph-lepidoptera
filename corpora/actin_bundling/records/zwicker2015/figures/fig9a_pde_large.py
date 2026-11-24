"""Calibrated reader for Fig. 9(a) autocatalytic PDE coarsening rates for large system (orange squares, V_s/w^3 = 8.8e5)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig9


def read(image, at):
    return _fig9.read_fig9a(image, at, "pde_large")
