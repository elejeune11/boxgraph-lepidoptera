"""Calibrated reader for Fig. 3A power-law fit line L ~ t^alpha."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3a import read_length_fit


def read(image, at):
    return read_length_fit(image, at)
