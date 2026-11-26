"""Calibrated reader for Fig. 5D (p. 8): barbed-end dwell time tau_1 for hProI^H133S (1, 2 uM AP)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig5de_bars(image, [0.0, 0.5, 1.0, 1.5, 2.0], 'hpro1_H133S')
