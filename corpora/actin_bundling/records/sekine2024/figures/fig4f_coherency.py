"""Digitize 4-M coherency (blue circles, 0-0.3) vs time (0-40 s) in Fig. 4f (PDF p. 7, 1200 dpi)."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4f_common import read_fig4f


def read(image, at):
    return read_fig4f(image, "coherency")
