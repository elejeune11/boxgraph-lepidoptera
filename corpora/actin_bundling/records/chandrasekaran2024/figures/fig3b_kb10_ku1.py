"""Thin wrapper for fig3b_kb10_ku1."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import _fig3b


def read(image, at):
    return _fig3b.read_fig3b(image, 0, '1.0', at)
