"""Calibrated reader for Fig. 6B (p. 9): subunits added before capping vs [CP] with 50 nM VASP (open circles)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_fig6b(image)
    pts = C.detect_open_circles(
        image, xa, ya, fr,
        exclude_fn=lambda x, y: (x > 6.0 and y > 4800.0) or (x > 2.5 and y > 5500.0) or x < 0.0,
        thr=228.0, close_k=7, area_range=(150, 620), dim_range=(14, 34),
    )
    nominal = [0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 20.0]
    if len(pts) != len(nominal):
        raise ValueError(f"Fig. 6B +VASP: expected {len(nominal)} open circles, found {len(pts)}: {[p['x'] for p in pts]}")
    for cp, p in zip(nominal, pts):
        p["x"] = cp
    return pts
