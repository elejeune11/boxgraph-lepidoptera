"""Calibrated figure reader for fig2a_rhoB_static."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig1_3 import read_fig2a

def read(image, at):
    return read_fig2a(image, at, "rhoB_static")
