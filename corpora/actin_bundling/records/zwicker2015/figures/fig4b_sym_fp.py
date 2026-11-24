"""Calibrated reader for Fig. 4(b) active symmetric stable two-droplet fixed point (R_sym/w, R_sym/w)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig4


def read(image, at):
    return _fig4.read_fig4_sym_fp(image, filled=True)
