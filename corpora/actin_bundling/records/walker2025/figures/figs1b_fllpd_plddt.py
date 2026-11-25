"""Calibrated reader for figs1b_fllpd_plddt: Fig. S1B full-length Lamellipodin (FL Lpd) SPARROW predicted pLDDT score."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _sparrow as SP

def read(image, at):
    return SP.read_sparrow_panel(image, "plddt", [200, 400, 600, 800, 1000, 1200], [float(v) for v in at], label="Fig. S1B full-length Lamellipodin (FL Lpd)")
