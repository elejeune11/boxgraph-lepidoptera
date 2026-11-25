"""Calibrated reader for Fig. 1E actin (green) transverse fluorescence intensity profile."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig1e import read_curve


def read(image, at):
    return read_curve(image, at, "actin")
