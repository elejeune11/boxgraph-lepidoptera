"""Fig. S1 Active experiments: growth exponent alpha vs. Lo domain area fraction phi_A."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS1 import read_series


def read(image, at):
    return read_series(image, "active")
