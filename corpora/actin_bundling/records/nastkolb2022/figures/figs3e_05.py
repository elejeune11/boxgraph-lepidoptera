"""Digitize Fig. S3E raw TIRFM actin kinetics at 0.5 uM actin, 0.25 uM VASP (dark gray squares; PDF p. 17, 1200 dpi)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs3e_common import read_figs3e

def read(image, at):
    return read_figs3e(image, "05", at)
