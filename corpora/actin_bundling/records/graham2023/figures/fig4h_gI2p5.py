"""Calibrated reader for Main Fig. 4h: gamma_I = 2.5 pN critical ring thickness."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4h import read_gi


def read(image, at):
    return read_gi(image, at, "gI2p5")
