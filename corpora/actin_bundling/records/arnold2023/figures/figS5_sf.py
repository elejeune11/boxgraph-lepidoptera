"""Fig. S5 Structure factor (black) curve: passive Cahn-Hilliard a/a_0 vs. t."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS5 import read_series


def read(image, at):
    return read_series(image, at, "sf")
