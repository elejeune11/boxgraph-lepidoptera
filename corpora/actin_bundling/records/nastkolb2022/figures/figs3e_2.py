"""Digitize Fig. S3E raw TIRFM actin kinetics at 2 uM actin, 0.25 uM VASP (blue up-triangles; PDF p. 17, 1200 dpi)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs3e_common import read_figs3e

def read(image, at):
    return read_figs3e(image, "2", at)
