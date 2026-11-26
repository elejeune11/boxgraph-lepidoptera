"""Calibrated reader for Fig. 5C — AFINES short cross-linker domain length vs k_grow."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _panels


def read(image, at):
    return _panels.read_fig5c(image, at, "short")
