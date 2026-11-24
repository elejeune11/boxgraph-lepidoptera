"""Calibrated reader for SI Fig. S6(B) right (High actin k_m^2 S(k,t)), t100."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS6 import read_s6_curve

def read(image, at):
    return read_s6_curve(image, at, "b_right", "t100")
