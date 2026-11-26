"""Calibrated reader for Fig. 3C Left — experimental fascin domain length (Ref. 16) vs [fascin]/[alpha-actinin]."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _panels


def read(image, at):
    return _panels.read_fig3c_left_expt(image, at)
