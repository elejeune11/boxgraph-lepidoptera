"""Calibrated reader for Supplementary Fig. 7b: Actin fluorescence intensity (a.u.) line profile across a single bulk fascin-actin bundle at fascin/actin = 0.5 vs Distance (um) (black curve illustrating width Delta)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_sfig7b(image, at)
