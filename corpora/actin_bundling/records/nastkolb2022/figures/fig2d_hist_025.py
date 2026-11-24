"""Digitize Fig. 2D bundle width histogram markers at 0.25 uM VASP (blue up-triangles; PDF p. 4, 1200 dpi)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2d_common import read_fig2d_hist

def read(image, at):
    return read_fig2d_hist(image, "025", at)
