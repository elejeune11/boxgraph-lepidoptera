"""Calibrated reader for Fig. 7(b) infinite-system critical radius R_crit/w vs phibar (dotted orange curve)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig7


def read(image, at):
    return _fig7.read_fig7b(image, at, "inf_Rcrit")
