"""Calibrated reader for Fig. 4(a) passive two-droplet blue nullcline (dR_1/dt = 0) R_2/w vs R_1/w."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig4


def read(image, at):
    return _fig4.read_fig4_nullcline_R1(image, at)
