"""Calibrated reader for Main Fig. 6f (Droplet 2 VASP normalized total intensity)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig6f import read_series


def read(image, at):
    return read_series(image, at, "d2_vasp")
