"""Calibrated reader for Fig. 1 inset: asymptotic lambda_tot (horizontal dashed line)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _henle_figs import read_fig1_inset


def read(image, at):
    return read_fig1_inset(image, at, "lam_tot_asymp")
