"""Calibrated figure reader for figS3_t60_k0007 (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_figS3_hist(image, at, x_bounds=(2.0, 509.0), y_bounds=(0.0, 100.0), rgb=C.MATLAB_ORANGE, x_scale="log10", x_interior=[4,9,17,34,66,131,258], y_interior=[20,40,60,80])
