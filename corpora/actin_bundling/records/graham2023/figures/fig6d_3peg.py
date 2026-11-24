"""Calibrated reader for Main Fig. 6d (3% PEG morphology abundance)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _stacked_bars import read_stacked_bars


def read(image, at):
    return read_stacked_bars(image, n_boxes=2, box_idx=0, n_bars=4, label="Fig6d 3%PEG")
