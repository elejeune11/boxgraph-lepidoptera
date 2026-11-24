"""Digitize Fig. 3E normalized TIRFM actin kinetics for 0.25 uM VASP +Profilin,CP (blue up-triangles; PDF p. 5, 1200 dpi)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3e_common import read_fig3e

def read(image, at):
    return read_fig3e(image, "vasp", at)
