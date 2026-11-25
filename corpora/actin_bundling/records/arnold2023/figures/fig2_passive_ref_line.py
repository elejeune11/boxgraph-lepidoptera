"""Fig. 2 Passive scaling reference line (black dashed line, alpha = 1/3)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2 import read_series


def read(image, at):
    return read_series(image, at, "passive_ref_line")
