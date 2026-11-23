"""Digitize Stage 15 PIV longitudinal mean displacement Vx (red, mean +/- SD) in Fig. 5c (PDF p. 9, 1200 dpi)."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5c_common import read_fig5c


def read(image, at):
    return read_fig5c(image, "vx")
