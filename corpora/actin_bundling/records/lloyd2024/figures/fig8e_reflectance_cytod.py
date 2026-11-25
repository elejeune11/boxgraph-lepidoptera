"""Digitize Fig. 8E: cytochalasin D treated average reflectance spectrum +- SE (blue), p. 11, 1200 dpi."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig8 import read_fig8e_spectrum


def read(image, at):
    return read_fig8e_spectrum(image, at, series="cytod")
