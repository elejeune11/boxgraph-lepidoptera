"""Calibrated reader for Fig. 5E Bottom — AFINES long cross-linker domain length vs k_grow (k_short^off,1/k_long^off,1 = 4)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _panels


def read(image, at):
    return _panels.read_fig5e(image, at, "bot", "long")
