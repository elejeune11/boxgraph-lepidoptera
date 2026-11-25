"""Digitize SI Fig. 1G: iridescent blue H. sara crossrib correlation function G_x(z) (p. 23, 1200 dpi)."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs1 import read_corr_panel


def read(image, at):
    # Printed y-ticks: 1.000, 0.999, 0.999 (two identical 3-dp rounded labels).
    # Equal spacing implies y_bot in [0.9985, 0.9990) so midpoint (1+y_bot)/2 in [0.99925, 0.9995) rounds to 0.999;
    # use midpoint estimate y_bot = 0.99875 with exact half-width uncertainty +-0.00025.
    return read_corr_panel(image, at, y_top=1.000, y_bot=0.99875, y_bot_uncert=0.00025)
