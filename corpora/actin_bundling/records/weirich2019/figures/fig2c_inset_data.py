"""Calibrated reader for Fig. 2C Inset — single myosin punctum trajectory dots d(t)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2c_inset import read_data


def read(image, at):
    return read_data(image)
