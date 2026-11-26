"""Calibrated reader for Fig. S3A Left — KMC waiting-time distribution before binding same cross-linker."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _panels


def read(image, at):
    return _panels.read_figS3a_left(image, at, "same")
