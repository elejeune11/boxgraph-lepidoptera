"""Calibrated reader for Fig. 5E (p. 8): actin monomers delivered per VASP tetramer for hProI^H133S (1, 2 uM AP)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig5de_bars(image, [0.0, 10.0, 20.0, 30.0, 40.0], 'hpro1_H133S')
