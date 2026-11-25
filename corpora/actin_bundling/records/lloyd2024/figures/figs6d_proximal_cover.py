"""Digitize SI Fig. 6D: 5 proximal cover scale markers (blue circles), p. 29, 1200 dpi."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs6 import read_figs6d_markers


def read(image, at):
    return read_figs6d_markers(image, region="proximal", scale_type="cover")
