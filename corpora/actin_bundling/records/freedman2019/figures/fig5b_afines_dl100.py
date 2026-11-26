"""Calibrated reader for Fig. 5B — AFINES short cross-linker domain length (Delta l_xl = 100 nm) vs rho_short/rho_long."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _panels


def read(image, at):
    return _panels.read_fig5b(image, at, "dl100")
