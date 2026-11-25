"""Fig. 3(d) Shear flow: Cahn-Hilliard domain size a(t) vs. t (open circles + fit line, alpha = 0.73)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3d import read_series


def read(image, at):
    return read_series(image, at, "shear")
