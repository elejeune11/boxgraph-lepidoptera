"""Calibrated reader for Fig. 9(b) red dashed horizontal line Q_stab/(w^3 k_0) from Eq. (25) with phi_0^B = 0."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig9


def read(image, at):
    return _fig9.read_fig9b(image, at, "Qstab_line")
