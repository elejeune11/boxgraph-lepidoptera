"""Calibrated reader for Fig. 5-S1 III(e): athermal CytoSim relaxation curve S(t) for S_0 = 0.1."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5s1e import read_sim


def read(image, at):
    return read_sim(image, at, 0.1)
