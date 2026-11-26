"""Calibrated reader for Fig. 1 inset: lambda_tot(R) for n_q = 10 uM (solid curve)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _henle_figs import read_fig1_inset


def read(image, at):
    return read_fig1_inset(image, at, "lam_tot")
