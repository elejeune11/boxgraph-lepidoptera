"""Calibrated reader for Supplementary Fig. 7d: Mean actin bundle width Delta (um) of encapsulated central aster arms at fascin/actin = 0.1, alpha-actinin/actin = 0.3 (red bar + SD error bar, 20 bundles in 10 GUVs)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_sfig7d(image, "aster")
