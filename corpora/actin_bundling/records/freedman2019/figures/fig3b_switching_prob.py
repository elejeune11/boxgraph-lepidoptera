"""Calibrated reader for Fig. 3B — conditional dissimilar cross-linker switching probability vs gap sites."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _panels


def read(image, at):
    return _panels.read_fig3b(image, at)
