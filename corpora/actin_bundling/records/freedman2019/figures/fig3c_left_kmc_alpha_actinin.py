"""Calibrated reader for Fig. 3C Left — KMC alpha-actinin domain length vs [fascin]/[alpha-actinin]."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _panels


def read(image, at):
    return _panels.read_fig3c_left_kmc(image, at, "alpha_actinin")
