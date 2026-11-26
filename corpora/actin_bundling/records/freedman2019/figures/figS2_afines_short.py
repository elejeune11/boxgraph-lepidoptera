"""Calibrated reader for Fig. S2 — AFINES calibration short cross-linker domain length vs shifted Delta l_xl."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _panels


def read(image, at):
    return _panels.read_figS2(image, at, "short")
