"""Fig. S6(A): simulated phase-field mass <phi>(t) for experimental PIV flow (Pe = 10^-4)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS6 import read_panel


def read(image, at):
    return read_panel(image, at, [0.0, 25.0, 50.0, 75.0, 100.0], "A")
