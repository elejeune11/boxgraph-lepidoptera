"""Calibrated reader for Fig. 3b: Actin fluorescence intensity (a.u.) line scan across GUV vs Length/GUV diameter at alpha-actinin/actin = 0.3, fascin/actin = 0.1 (black curve)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig3b(image, at, "0.3")
