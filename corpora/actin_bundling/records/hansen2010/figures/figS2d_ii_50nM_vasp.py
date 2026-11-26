"""Calibrated reader for Fig. S2D(ii) (p. 16): K_obs vs [CP] at 100 mM KCl with 50 nM VASP (open circles)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_figS2d_ii(image)
    pts = C.detect_open_circles(
        image, xa, ya, fr,
        exclude_fn=lambda x, y: (x > 4.5 and y > 0.038) or (x > 3.2 and y < 0.011),
        thr=230.0, close_k=7, area_range=(150, 500), dim_range=(12, 28),
    )
    if len(pts) != 4:
        raise ValueError(f"Fig. S2D(ii) +VASP: expected 4 open circles, found {len(pts)}: {[p['x'] for p in pts]}")
    return pts
