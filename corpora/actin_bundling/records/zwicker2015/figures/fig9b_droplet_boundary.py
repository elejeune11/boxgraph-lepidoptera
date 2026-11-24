"""Calibrated reader for Fig. 9(b) black boundary separating no-droplets from droplet-forming regions."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig9


def read(image, at):
    return _fig9.read_fig9b(image, at, "droplet_boundary")
