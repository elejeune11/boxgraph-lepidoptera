"""Calibrated reader for Fig. 4(C) left (Experiment): No Actin curvature distribution P(kappa)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4 import read_4c_curve

def read(image, at):
    return read_4c_curve(image, at, "exp", "no_actin")
