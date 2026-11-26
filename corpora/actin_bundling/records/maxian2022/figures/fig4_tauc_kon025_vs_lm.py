"""Calibrated figure reader for fig4_tauc_kon025_vs_lm (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig4_heatmap_col(image, 0.25)
