"""Calibrated reader for Fig. 2B — myosin probability distribution at t > 5 min."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2b import read_series


def read(image, at):
    return read_series(image, "t_gt_5min")
