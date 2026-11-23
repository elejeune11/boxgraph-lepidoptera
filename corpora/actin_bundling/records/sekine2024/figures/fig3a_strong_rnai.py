"""Digitize 10 blue strong-defect RNAi markers in Fig. 3a (PDF p. 5, 900 dpi)."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3a_common import read_fig3a_markers


def read(image, at):
    return read_fig3a_markers(image, "blue")
