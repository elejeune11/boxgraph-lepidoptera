"""Calibrated reader for Fig. 3A constant aspect-ratio fit line."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3a import read_aspect_ratio_fit


def read(image, at):
    return read_aspect_ratio_fit(image, at)
