"""Calibrated reader for Fig. S3B (p. 27): 647-BSA flow-in RAM vs time."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_figS3b(image)
