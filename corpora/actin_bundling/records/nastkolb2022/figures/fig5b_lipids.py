"""Digitize Fig. 5B normalized lipid fluorescence cross-section (blue curve, left axis; PDF p. 8, 1200 dpi)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5b_common import read_fig5b

def read(image, at):
    return read_fig5b(image, "lipids", at)
