"""Calibrated reader for Fig. S3 black filled circle series (R_i = 2.2 um)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS3 import read_series


def read(image, at):
    return read_series(image, "Ri_2p2")
