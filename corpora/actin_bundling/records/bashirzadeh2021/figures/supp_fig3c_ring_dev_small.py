"""Calibrated reader for Supplementary Fig. 3c: Normalized radial deviation d/d_avg vs normalized arc length s/p for an actin ring in a small GUV (blue curve)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_sfig3c(image, at, "small")
