"""Calibrated figure reader for figS4h_freq_sim."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS4_S5 import read_figS4h_sim

def read(image, at):
    return read_figS4h_sim(image, at)
