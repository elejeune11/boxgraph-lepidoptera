"""Fig. S6(B): simulated phase-field mass <phi>(t) for rotational flow."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS6 import read_panel


def read(image, at):
    return read_panel(image, at, [0.0, 20.0, 40.0], "B")
