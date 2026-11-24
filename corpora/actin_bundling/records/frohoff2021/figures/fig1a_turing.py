"""Calibrated reader for Fig. 1a (turing boundary)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig1 import read_fig1_curve

def read(image, at):
    return read_fig1_curve(image, at, '1a', 'turing')
