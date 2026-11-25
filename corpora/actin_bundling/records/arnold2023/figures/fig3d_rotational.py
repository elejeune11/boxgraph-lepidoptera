"""Fig. 3(d) Rotational flow: Cahn-Hilliard domain size a(t) vs. t (open squares + fit line, alpha = 0.28)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3d import read_series


def read(image, at):
    return read_series(image, at, "rotational")
