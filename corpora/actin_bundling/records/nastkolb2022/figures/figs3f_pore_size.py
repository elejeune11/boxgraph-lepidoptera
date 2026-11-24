"""Digitize Fig. S3F pore size boxplots (median + 1-99th whiskers) vs actin concentration (PDF p. 17, 1200 dpi)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs3f_common import read_figs3f

def read(image, at):
    return read_figs3f(image, "pore")
