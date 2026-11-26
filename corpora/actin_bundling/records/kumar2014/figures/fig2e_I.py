"""Calibrated reader for Fig. 2(e) normalized inhibitor profile I/I_0 at time t+T/2."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2 import read_fig2_profile


def read(image, at):
    return read_fig2_profile(image, at, "I", "e")
