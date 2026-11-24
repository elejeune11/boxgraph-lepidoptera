"""Calibrated reader for SI Fig. S6(A) left (Low actin S(k,t)), t100."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS6 import read_s6_curve

def read(image, at):
    return read_s6_curve(image, at, "a_left", "t100")
