"""Calibrated reader for Fig. 3 (top): No actin (alpha = 1/3 dashed reference line)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3 import read_sim_no_actin

def read(image, at):
    return read_sim_no_actin(image, at)
