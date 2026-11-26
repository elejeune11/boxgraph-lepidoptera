"""Calibrated reader for Fig. 4F (p. 5): filament assembly boundary line RAM*c = 0.1 uM."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig4f(image, at, "filament_boundary")
