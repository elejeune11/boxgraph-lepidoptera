"""Calibrated reader for figs3d_frap_actin: Fig. S3D actin FRAP (+/- SD)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _frap_si as FS

def read(image, at):
    return FS.read_frap_si(image, "actin", [0, 50, 100, 150, 200, 250], [float(v) for v in at], label="Fig. S3D")
