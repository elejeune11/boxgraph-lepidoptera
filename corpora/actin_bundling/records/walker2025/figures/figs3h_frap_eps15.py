"""Calibrated reader for figs3h_frap_eps15: Fig. S3H Eps15 FRAP (+/- SD)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _frap_si as FS

def read(image, at):
    return FS.read_frap_si(image, "protein", [0, 50, 100, 150, 200], [float(v) for v in at], label="Fig. S3H")
