"""Fig. S4 Early time (purple): near-critical passive bilayer a/a_0 vs. t [s]."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS4 import read_series


def read(image, at):
    return read_series(image, at, "early")
