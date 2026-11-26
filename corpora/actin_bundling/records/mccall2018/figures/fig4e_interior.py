"""Calibrated reader for Fig. 4E (p. 5): mean droplet interior intensity vs [TMR-actin]."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig4e(image, "interior")
