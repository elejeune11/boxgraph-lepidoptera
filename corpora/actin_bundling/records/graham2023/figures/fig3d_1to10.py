"""Calibrated reader for Main Fig. 3d (1:10 shell/ring/disc abundance)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _stacked_bars import read_stacked_bars


def read(image, at):
    return read_stacked_bars(image, n_boxes=3, box_idx=1, n_bars=3, label="Fig3d 1:10")
