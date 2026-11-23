"""Calibrated reader for Fig. 3A Inset (p. 3): magnified SAXS diffraction I(q) at R* = 0.5."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3a import read_fig3a_series


def read(image, at):
    return read_fig3a_series(image, at, "inset_R0p5")
