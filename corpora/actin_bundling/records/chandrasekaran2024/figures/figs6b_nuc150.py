"""Thin wrapper for figs6b_nuc150."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import _figs6b


def read(image, at):
    return _figs6b.read_figs6b(image, 7, at)
