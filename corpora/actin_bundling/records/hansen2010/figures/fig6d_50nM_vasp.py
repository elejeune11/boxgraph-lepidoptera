"""Calibrated reader for Fig. 6D (p. 9): K_obs vs [CP] with 50 nM VASP (open white circles)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_fig6d(image)
    pts = C.detect_open_circles(
        image, xa, ya, fr,
        exclude_fn=lambda x, y: (x > 7.5 and y > 0.036) or (x > 7.5 and y < 0.012),
        thr=236.0, close_k=7, area_range=(200, 620), dim_range=(14, 34),
    )
    if len(pts) != 5:
        raise ValueError(f"Fig. 6D +VASP: expected 5 open circles, found {len(pts)}: {[p['x'] for p in pts]}")
    return pts
