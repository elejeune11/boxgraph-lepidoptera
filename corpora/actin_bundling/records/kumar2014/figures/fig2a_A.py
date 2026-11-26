"""Calibrated reader for Fig. 2(a) normalized activator profile A/A_0 at time t."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2 import read_fig2_profile


def read(image, at):
    return read_fig2_profile(image, at, "A", "a")
