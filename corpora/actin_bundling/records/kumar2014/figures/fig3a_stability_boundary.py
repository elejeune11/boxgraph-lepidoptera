"""Calibrated reader for Fig. 3(a) 2D ARD homogeneous-to-patterned stability onset P_c(beta)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3a import read_stability_boundary


def read(image, at):
    return read_stability_boundary(image, at)
