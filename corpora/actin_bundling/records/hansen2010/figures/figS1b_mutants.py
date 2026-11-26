"""Calibrated reader for Fig. S1B (p. 15): horizontal bar chart of polymerization rate for GAB/FAB mutants."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_figS1b_bars(image)
