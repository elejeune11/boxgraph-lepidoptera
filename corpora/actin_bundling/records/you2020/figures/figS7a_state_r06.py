"""Calibrated figure reader for figS7a_state_r06."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS7 import read_figS7a_state_r06

def read(image, at):
    return read_figS7a_state_r06(image, at)
