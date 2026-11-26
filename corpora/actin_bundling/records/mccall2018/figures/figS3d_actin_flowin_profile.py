"""Calibrated reader for Fig. S3D (p. 27): TMR-actin intensity line scan 15 s post flow-in."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_figS3d(image, at, "profile")
