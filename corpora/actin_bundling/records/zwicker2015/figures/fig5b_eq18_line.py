"""Calibrated reader for Fig. 5(b) analytical passive coarsening rate Eq. (18) (legend prints Eq. 16) vs V_bar/w^3."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig5


def read(image, at):
    return _fig5.read_fig5b(image, at, "line")
