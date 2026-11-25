"""Digitize SI Fig. 4E: per-ridge mean width +- SD across Ridges 1..30 (PDF p. 26, 1200 dpi)."""

from __future__ import annotations

from pathlib import Path
import sys
import numpy as np
from boxgraph import digitize as D

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import run_centers, y_interval


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    gray = image.astype(float).mean(axis=2)
    H, W = gray.shape

    c_cands = np.where((gray[100:2650, 200:800] < 215.0).mean(axis=0) > 0.85)[0] + 200
    c0 = float(c_cands.mean() + 0.5)
    r_cands = np.where((gray[2500:2850, 600:4400] < 215.0).mean(axis=1) > 0.85)[0] + 2500
    r1 = float(r_cands.mean() + 0.5)

    yt = run_centers(np.where(gray[150:int(r1) - 10, int(round(c0)) - 8] < 210.0)[0] + 150, gap=4)
    xt = run_centers(np.where(gray[int(round(r1)) + 8, int(round(c0)) + 10:W - 50] < 210.0)[0] + int(round(c0)) + 10, gap=4)
    if len(yt) != 5 or len(xt) != 30:
        raise RuntimeError(f"SI Fig. 4E expected 5 y-ticks and 30 x-ticks, got yt={len(yt)}, xt={len(xt)}")

    ya = D.Axis("linear", list(zip(yt, [0.5, 0.4, 0.3, 0.2, 0.1])))
    # Jet-black mean dots and vertical SD lines have gray < 65 (grey individual points are > 110)
    blk = gray < 65.0

    out = []
    for rid, xc in enumerate(xt, start=1):
        c_lo, c_hi = int(round(xc)) - 30, int(round(xc)) + 31
        sub = blk[100:int(r1) - 20, c_lo:c_hi]
        rows_b = np.where(sub.any(axis=1))[0] + 100
        widths = sub.sum(axis=1)
        dot_rows = np.where(widths >= 20)[0] + 100
        if len(dot_rows) == 0 or len(rows_b) == 0:
            out.append({"x": float(rid), "y": None, "note": f"ridge {rid} black dot/SD not found"})
            continue
        r_mean = float(dot_rows.mean() + 0.5)
        r_top = float(rows_b.min() + 0.5)
        r_bot = float(rows_b.max() + 0.5)

        y_mean, lo, hi, _, _ = y_interval(ya, r_mean, 2.5)
        bar_hi = float(ya.to_value(r_top))
        bar_lo = float(ya.to_value(r_bot))
        h_px = 3.0 + max(ya.residual_px, 0.5)
        bar_hi_outer = float(ya.to_value(r_top - h_px))
        bar_lo_outer = float(ya.to_value(r_bot + h_px))
        out.append({
            "x": float(rid),
            "y": float(y_mean), "lo": float(lo), "hi": float(hi),
            "bar_lo": float(min(bar_lo, bar_hi)), "bar_hi": float(max(bar_lo, bar_hi)),
            "bar_lo_outer": float(min(bar_lo_outer, bar_hi_outer)),
            "bar_hi_outer": float(max(bar_lo_outer, bar_hi_outer)),
            "px": [float(xc), r_mean],
            "px_lo": float(r_top), "px_hi": float(r_bot),
            "note": f"Ridge {rid} mean +- SD (ry={ya.residual_px:.2f} px)",
        })
    return out
