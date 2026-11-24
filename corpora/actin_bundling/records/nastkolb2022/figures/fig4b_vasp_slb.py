"""Digitize Fig. 4B Spearman colocalization score vs time for 0.25 uM VASP on SLB (black squares; PDF p. 7, 1200 dpi)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4b_common import read_fig4b

def read(image, at):
    return read_fig4b(image, "vasp_slb")
