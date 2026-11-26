"""Calibrated reader for Fig. S4C — TIRF measured fascin domain length with unlabeled actin."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _panels


def read(image, at):
    return _panels.read_figS4c(image, at, "fascin")
