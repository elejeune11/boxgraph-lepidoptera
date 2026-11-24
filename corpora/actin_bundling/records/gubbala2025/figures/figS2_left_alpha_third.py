"""Calibrated reader for SI Fig. S2 (left): alpha = 1/3 dashed reference line."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS2 import read_left_alpha_third

def read(image, at):
    return read_left_alpha_third(image, at)
