"""Calibrated figure reader for fig2b_v_sim."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig1_3 import read_fig2b

def read(image, at):
    return read_fig2b(image, at, "sim")
