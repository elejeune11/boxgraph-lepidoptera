"""Calibrated reader for Fig. 3 (top): High actin density simulation (light blue/teal circles & line)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3 import read_sim_colored

def read(image, at):
    return read_sim_colored(image, at, "high")
