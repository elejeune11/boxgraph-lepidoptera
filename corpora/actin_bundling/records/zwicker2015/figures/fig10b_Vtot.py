"""Calibrated reader for Fig. 10(b) total steady-state droplet volume (V_1 + V_2)/w^3 vs Q_1/(w^3 k_0) (dotted black curve)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig10


def read(image, at):
    return _fig10.read_fig10b(image, at, "Vtot")
