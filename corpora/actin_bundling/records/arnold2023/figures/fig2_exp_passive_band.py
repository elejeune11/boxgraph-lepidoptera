"""Fig. 2 Passive experiments (combined band of 3 overlapping black closed-symbol datasets)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2 import read_series


def read(image, at):
    return read_series(image, at, "exp_passive_band")
