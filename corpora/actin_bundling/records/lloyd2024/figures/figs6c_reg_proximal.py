"""Digitize SI Fig. 6C: H. e. demophoon proximal (blue) regression +- 95% CI (ridge number vs scale width), p. 29, 1200 dpi."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs6 import read_figs6c_reg


def read(image, at):
    return read_figs6c_reg(image, at, series="blue")
