"""Calibrated reader for Fig. 1(e) numerical PDE oscillation period markers T/tau vs beta."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig1e import read_pde_period


def read(image, at):
    return read_pde_period(image, at)
