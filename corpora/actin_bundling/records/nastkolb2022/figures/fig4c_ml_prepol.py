"""Digitize Fig. 4C FRAP recovery curve for VASP ML prePol (blue; PDF p. 7, 1200 dpi)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4c_common import read_fig4c

def read(image, at):
    return read_fig4c(image, "ml_prepol", at)
