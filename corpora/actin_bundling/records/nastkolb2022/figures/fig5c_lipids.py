"""Digitize Fig. 5C Lipids (Texas Red-DHPE) intensity ratio under/no-bundles (blue curve +/- SD, n=5 bundles; PDF p. 8, 1200 dpi)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5c_common import read_fig5c

def read(image, at):
    return read_fig5c(image, "lipids", at)
