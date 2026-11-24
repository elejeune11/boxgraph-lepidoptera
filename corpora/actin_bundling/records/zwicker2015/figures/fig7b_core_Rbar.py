"""Calibrated reader for Fig. 7(b) finite-system stable radius R_bar/w vs phibar with active core Q_1 = 300 w^3 k_0 (solid green curve)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig7


def read(image, at):
    return _fig7.read_fig7b(image, at, "core_Rbar")
