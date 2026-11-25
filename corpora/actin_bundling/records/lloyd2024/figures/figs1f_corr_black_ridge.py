"""Digitize SI Fig. 1F: non-iridescent black H. sara ridge correlation function G_z(y) (p. 23, 1200 dpi)."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs1 import read_corr_panel


def read(image, at):
    # Printed y-ticks: 1.000, 0.992, 0.984 (3-dp rounded linspace -> bottom tick 0.984 +- 0.0005)
    return read_corr_panel(image, at, y_top=1.000, y_bot=0.984, y_bot_uncert=0.0005)
