"""Calibrated reader for Extended Data Fig. 4c (1:7.5 no-EDTA/EGTA shell/ring/disc)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _stacked_bars import read_stacked_bars


def read(image, at):
    return read_stacked_bars(image, n_boxes=3, box_idx=2, n_bars=3, label="ED4c 1:7.5")
