"""Fig. 2 PIV Cahn-Hilliard theory Pe = 10^-4 (green + + dashed line, alpha = 0.34)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2 import read_series


def read(image, at):
    return read_series(image, at, "sim_piv_pe1e4")
