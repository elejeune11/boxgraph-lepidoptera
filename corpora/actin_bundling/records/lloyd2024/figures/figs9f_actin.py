"""Calibrated reader for SI Fig. 9F green actin (phalloidin) intensity profile (PDF p. 32, 1200 dpi)."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fluo_profiles import read_figs9f


def read(image, at=None, **kwargs):
    xs = [float(v) for v in at] if at is not None else [round(0.15 * i, 2) for i in range(87)]
    return read_figs9f(image, channel="actin", at=xs)
