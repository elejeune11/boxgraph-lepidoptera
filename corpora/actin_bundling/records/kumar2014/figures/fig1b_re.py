"""Calibrated reader for Fig. 1(b) Re[tau*lambda_+] vs k*ell."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig1abc import read_dispersion


def read(image, at):
    return read_dispersion(image, at, "re", "b")
