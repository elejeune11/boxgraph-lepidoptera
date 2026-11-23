"""Digitize simulated S(q) curve (blue markers) in Supplementary Fig. 7c4 (a4-1 low-turnover clusters, PDF p. 22, 1200 dpi)."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _s7c_common import read_s7c_sim


def read(image, at):
    return read_s7c_sim(image, at, [0.0, 2.0, 4.0, 6.0], [3.0, 2.0, 1.0, 0.0])
