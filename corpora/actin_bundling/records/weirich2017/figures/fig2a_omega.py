"""Calibrated reader for Fig. 2A anchoring strength omega = gamma_A/gamma_I vs filamin %."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2a import read_series


def read(image, at):
    return read_series(image, "omega")
