"""Calibrated reader for Fig. S2D(i) (p. 16): subunits added before capping at 100 mM KCl with 50 nM VASP (open circles)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_figS2d_i(image)
    pts = C.detect_open_circles(
        image, xa, ya, fr,
        exclude_fn=lambda x, y: (x > 5.0 and y > 4000.0) or (x > 1.0 and y > 5000.0),
        thr=228.0, close_k=7, area_range=(150, 500), dim_range=(12, 28),
    )
    nominal = [0.0, 1.0, 2.0, 4.0, 8.0, 12.0]
    if len(pts) != len(nominal):
        raise ValueError(f"Fig. S2D(i) +VASP: expected 6 open circles, found {len(pts)}: {[p['x'] for p in pts]}")
    for cp, p in zip(nominal, pts):
        p["x"] = cp
    return pts
