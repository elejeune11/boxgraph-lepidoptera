"""Calibrated reader for Fig. S3D (p. 17): barbed-end polymerization rate at 1.5 uM actin for 4 constructs."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    labels = ["no VASP", "+VASP(1-380) LIL-3A/RRRK-4A", "+VASP(280-380aa)", "+GST-VASP(280-380aa)"]
    return C.read_vertical_gray_bars(image, [0, 5, 10, 15, 20], 4, labels)
