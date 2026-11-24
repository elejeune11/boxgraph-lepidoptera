"""Calibrated reader for Fig. 9 blue n+ subcriticality boundary Delta_+(Sigma)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig8_10 import read_fig9_subcrit_plus

def read(image, at):
    return read_fig9_subcrit_plus(image, at)
