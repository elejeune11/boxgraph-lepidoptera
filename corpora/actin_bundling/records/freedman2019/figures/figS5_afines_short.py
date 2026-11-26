"""Calibrated reader for Fig. S5 — AFINES short cross-linker domain length vs k_grow (1500-1600 s window)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _panels


def read(image, at):
    return _panels.read_figS5(image, at, "short")
