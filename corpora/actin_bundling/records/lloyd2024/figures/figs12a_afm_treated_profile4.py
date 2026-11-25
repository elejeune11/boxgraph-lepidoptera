"""Digitize SI Fig. 12A: Cyto-D treated transverse AFM height profile (Profile 4, blue), p. 36, 1800 dpi."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs12 import read_afm_profile


def read(image, at):
    return read_afm_profile(image, at, xmax=3.0, color="blue")
