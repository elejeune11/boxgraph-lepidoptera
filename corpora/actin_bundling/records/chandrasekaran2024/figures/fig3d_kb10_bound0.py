"""Thin wrapper for fig3d_kb10_bound0."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import _fig3d


def read(image, at):
    return _fig3d.read_fig3d(image, 0, 0)
