"""Calibrated reader for Fig. S4D (p. 18): Cy3-VASP fold increase on single filaments vs parallel/anti-parallel bundles."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    labels = ["single filament (n=30)", "parallel bundle (n=16)", "anti-parallel bundle (n=10)"]
    return C.read_vertical_gray_bars(image, [0, 1, 2, 3, 4, 5, 6, 7, 8], 3, labels)
