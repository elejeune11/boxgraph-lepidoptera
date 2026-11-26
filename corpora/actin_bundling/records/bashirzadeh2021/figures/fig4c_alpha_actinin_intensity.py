"""Calibrated reader for Fig. 4c: Normalized alpha-actinin fluorescence intensity along dashed line across central-aster GUV vs Length (um) (blue curve, sharply peaked at GUV center ~25 um)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig4c(image, at, "alpha_actinin")
