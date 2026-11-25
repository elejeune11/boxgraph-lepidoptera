"""Digitize Fig. 2D: adult H. sara blue scales linear regression +- 95% CI (ridge number vs scale width), p. 4, 1200 dpi."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2 import calibrate_fig2_scatter, read_lm_and_ci


def read(image, at):
    xa, ya, c0, r1 = calibrate_fig2_scatter(image, [25.0, 30.0, 35.0, 40.0, 45.0, 50.0], [60.0, 50.0, 40.0, 30.0])
    return read_lm_and_ci(image, xa, ya, c0, r1, r_top=400, c_right=2100,
                          series="blue", at=at, x_fit_range=(26.0, 48.2))
