"""Calibrated reader for Fig. 7E green actin (phalloidin) intensity profile (PDF p. 10, 1200 dpi)."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fluo_profiles import read_fig7ef


def read(image, at=None, **kwargs):
    xs = [float(v) for v in at] if at is not None else [round(0.1 * i, 2) for i in range(69)]
    return read_fig7ef(image, panel="7E", channel="actin", at=xs)
