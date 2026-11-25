"""Calibrated reader for figs1b_minilpd_disorder: Fig. S1B Lpd aa850-1250 SPARROW disorder score."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _sparrow as SP

def read(image, at):
    return SP.read_sparrow_panel(image, "disorder", [50, 100, 150, 200, 250, 300, 350, 400], [float(v) for v in at], label="Fig. S1B Lpd aa850-1250")
