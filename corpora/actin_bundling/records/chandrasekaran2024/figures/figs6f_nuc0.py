"""Thin wrapper for figs6f_nuc0."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import _figs6f


def read(image, at):
    return _figs6f.read_figs6f(image, 'nuc0', at)
