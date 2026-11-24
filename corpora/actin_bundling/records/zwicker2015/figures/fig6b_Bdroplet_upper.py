"""Calibrated reader for Fig. 6(b) upper B-droplet boundary k_b/k_0 = 16*(phi_0^B)^3 vs k_f/k_0."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig6


def read(image, at):
    return _fig6.read_fig6b(image, at, "Bdroplet_upper")
