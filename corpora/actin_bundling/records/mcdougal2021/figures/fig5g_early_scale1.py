"""Calibrated reader for Fig. 5G (p. 5): early individual scale 1 (top cyan curve, rgb 0,179,179)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5g import read_teal_scale

def read(image, at):
    return read_teal_scale(image, (0, 179, 179), tol=14.0)
