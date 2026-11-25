"""Fig. 3(c) Shear-flow static structure factor <S(k)> at t = 30 (squares)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3c import read_series


def read(image, at):
    return read_series(image, at, "t30")
