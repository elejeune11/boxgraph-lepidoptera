"""Calibrated reader for Fig. S3C — KMC effective switching barrier vs growth-to-bundling rate ratio (mu_grow = 20 nm/s)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _panels


def read(image, at):
    return _panels.read_figS3c(image, at, "20nms")
