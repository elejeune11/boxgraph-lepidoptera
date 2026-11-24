"""Calibrated reader for Fig. 7(a) finite-system autocatalytic droplet growth rate without core (blue curve)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig7


def read(image, at):
    return _fig7.read_fig7a(image, at, "fin")
