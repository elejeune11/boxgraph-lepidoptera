"""Calibrated reader for SI Fig. 10A treated (21 teal TREATED_*) peak reflectance (PDF p. 34, 900 dpi)."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs10 import read_figs10_peaks


def read(image, at=None, **kwargs):
    return read_figs10_peaks(image, panel="10A", group="treated")
