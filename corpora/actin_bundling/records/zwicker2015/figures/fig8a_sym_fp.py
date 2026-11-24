"""Calibrated reader for Fig. 8(a) autocatalytic symmetric unstable fixed point (Q_i = 0)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig8


def read(image, at):
    return _fig8.read_fig8_sym_fp(image, filled=False)
