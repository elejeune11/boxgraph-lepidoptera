"""Calibrated reader for Fig. 3D (p. 5): Lateral filament side dwell time tau_1 for +Actin (1) and +Actin, CP (2) with SEM error bars."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig3d_group(image, 'lateral_side')
