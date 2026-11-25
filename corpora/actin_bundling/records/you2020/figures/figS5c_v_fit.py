"""Calibrated figure reader for figS5c_v_fit."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS4_S5 import read_figS5c_fit

def read(image, at):
    return read_figS5c_fit(image, at)
