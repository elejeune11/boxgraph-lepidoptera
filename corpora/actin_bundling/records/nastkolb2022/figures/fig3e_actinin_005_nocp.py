"""Digitize Fig. 3E normalized TIRFM actin kinetics for 0.05 uM alpha-actinin1b -Profilin,CP (cyan left-triangles; PDF p. 5, 1200 dpi)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3e_common import read_fig3e

def read(image, at):
    return read_fig3e(image, "actinin_005_nocp", at)
