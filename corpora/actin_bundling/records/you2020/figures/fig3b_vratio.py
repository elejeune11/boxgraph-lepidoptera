"""Calibrated figure reader for fig3b_vratio."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig1_3 import read_fig3b

def read(image, at):
    return read_fig3b(image, at)
