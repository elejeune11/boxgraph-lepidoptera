"""Digitize Fig. 2E pore size boxplots (broken right y-axis, median + 1-99th whiskers) vs VASP concentration (PDF p. 4, 1200 dpi)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2e_common import read_fig2e

def read(image, at):
    return read_fig2e(image, "pore", include_whiskers=True)
