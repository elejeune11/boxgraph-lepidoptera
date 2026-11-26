"""Calibrated reader for Fig. 5B — experimental fascin domain length (Ref. 16) vs concentration ratio."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _panels


def read(image, at):
    return _panels.read_fig5b(image, at, "expt")
