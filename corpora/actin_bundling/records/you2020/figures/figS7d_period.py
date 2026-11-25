"""Calibrated figure reader for figS7d_period."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS7 import read_figS7d_period

def read(image, at):
    return read_figS7d_period(image, at)
