"""Fig. 2 Passive Cahn-Hilliard theory (green closed circles)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2 import read_series


def read(image, at):
    return read_series(image, at, "sim_passive")
