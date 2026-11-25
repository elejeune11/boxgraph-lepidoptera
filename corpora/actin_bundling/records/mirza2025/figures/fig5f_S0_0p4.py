"""Calibrated reader for Fig. 5(f): ensemble-averaged S(t)-S_0 vs t/tau_A for S_0 = 0.4 (green)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5 import read_fig5f_curve


def read(image, at):
    return read_fig5f_curve(image, at, 'S0_0p4')
