"""Calibrated reader for Supplementary Fig. 8e: Normalized alpha-actinin fluorescence intensity across encapsulated actin ring vs Length (um) at alpha-actinin/actin=0.3, fascin/actin=0.1 (magenta curve)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_sfig8e(image, at, "alpha_actinin")
