"""Calibrated reader for Fig. 3 (bottom): High actin density experiment (light blue/teal squares)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3 import read_exp_colored

def read(image, at):
    return read_exp_colored(image, at, "high")
