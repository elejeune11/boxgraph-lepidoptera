"""Digitize SI Fig. 12C: AFM ridge height mean +- SD for control (x=1) and cyto-D (x=2), p. 36, 1200 dpi."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs12 import read_s12c_mean_sd


def read(image, at):
    return read_s12c_mean_sd(image)
