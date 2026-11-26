"""Calibrated reader for Fig. 4B — TIRF measured alpha-actinin domain length at high and low labeled actin."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _panels


def read(image, at):
    return _panels.read_fig4b(image, at, "alpha_actinin")
