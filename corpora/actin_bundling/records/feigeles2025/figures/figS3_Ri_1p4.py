"""Calibrated reader for Fig. S3 red X series (R_i = 1.4 um)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS3 import read_series


def read(image, at):
    return read_series(image, "Ri_1p4")
