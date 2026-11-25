"""Fig. 2 Active experiment (black open squares, alpha = 0.69)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2 import read_series


def read(image, at):
    return read_series(image, at, "exp_active_069")
