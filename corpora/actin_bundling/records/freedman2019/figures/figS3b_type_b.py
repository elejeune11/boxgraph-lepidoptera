"""Calibrated reader for Fig. S3B — KMC bundle extension rate probability density (Type B)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _panels


def read(image, at):
    return _panels.read_figS3b(image, at)
