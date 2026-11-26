"""Calibrated reader for Fig. 3c: Cumulative % Periphery cluster formation in GUVs >20 um (1=alpha-Actinin alone n=87, 2=alpha-Actinin+fascin n=80; grey bars)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig3c(image, "periphery")
