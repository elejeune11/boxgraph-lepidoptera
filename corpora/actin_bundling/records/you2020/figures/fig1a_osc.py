"""Calibrated figure reader for fig1a_osc."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig1_3 import read_fig1a

def read(image, at):
    return read_fig1a(image, at, "osc")
