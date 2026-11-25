"""Fig. 2 Active experiment (black open diamonds, alpha = 0.59)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2 import read_series


def read(image, at):
    return read_series(image, at, "exp_active_059")
