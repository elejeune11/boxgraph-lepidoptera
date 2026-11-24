"""Digitize Fig. S4C log-normal fit curve for 0.25 uM VASP noCP on ML (orange curve; PDF p. 18, 1200 dpi)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs4c_common import read_figs4c_fit

def read(image, at):
    return read_figs4c_fit(image, "nocp_ml", at)
