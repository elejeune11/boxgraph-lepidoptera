"""Calibrated reader for SI Fig. S1 (p. 6): 8 translated experimental points (red filled circles)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS1 import read_figS1_experiment


def read(image, at):
    return read_figS1_experiment(image)
