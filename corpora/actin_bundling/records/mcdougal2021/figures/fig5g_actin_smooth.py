"""Calibrated reader for Fig. 5G (p. 5): 8 smooth-stage actin spacing black crosses (~38-44% dev)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5g import read_black_smooth

def read(image, at):
    return read_black_smooth(image)
