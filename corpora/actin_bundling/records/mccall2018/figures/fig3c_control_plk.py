"""Calibrated reader for Fig. 3C (p. 4): 1.5 uM actin control baseline (cyan dashed)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig3c(image, at, "control_plk")
