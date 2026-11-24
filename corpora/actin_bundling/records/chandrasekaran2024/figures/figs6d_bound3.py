"""Thin wrapper for figs6d_bound3."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import _figs6d


def read(image, at):
    return _figs6d.read_figs6d(image, 3)
