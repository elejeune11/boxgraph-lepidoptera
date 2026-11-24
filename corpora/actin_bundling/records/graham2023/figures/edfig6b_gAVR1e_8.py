"""Calibrated reader for Extended Data Fig. 6B (final energy along gamma_AVR = 1e-8 pN/nm)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _edfig6 import read_ed6_row


def read(image, at):
    return read_ed6_row(image, at, panel="B", row_idx=0)
