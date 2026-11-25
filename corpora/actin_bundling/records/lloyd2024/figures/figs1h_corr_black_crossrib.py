"""Digitize SI Fig. 1H: non-iridescent black H. sara crossrib correlation function G_x(z) (p. 23, 1200 dpi)."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs1 import read_corr_panel


def read(image, at):
    # Printed y-ticks: 1.000, 0.999, 0.998 (3-dp rounded linspace -> bottom tick 0.998 +- 0.0005)
    return read_corr_panel(image, at, y_top=1.000, y_bot=0.998, y_bot_uncert=0.0005)
