"""Calibrated reader for Fig. 3 (bottom inset): No actin experimental squares."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3 import read_inset_exp

def read(image, at):
    return read_inset_exp(image, at)
