"""Calibrated reader for Main Fig. 3B right (p. 4): ring spacing bar chart (WT, Src42A, kkv) with broken y-axis and error bars."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import clusters, verify_axis, value_interval


def read(image, at):
    gray = image.astype(float).mean(axis=2)
    dark = gray < 120.0
    H, W = dark.shape

    # Horizontal x-axis line (exclude bottom panel border near H)
    h_lines = [h for h in clusters(np.where(dark[: H - 80].mean(axis=1) > 0.45)[0]) if h[0] > 1000]
    if not h_lines:
        raise RuntimeError("Fig. 3B bar: could not locate x-axis baseline")
    r_base, rb0, rb1 = h_lines[-1]

    # Vertical y-axis spine column (left third of crop)
    v_col = int(np.argmax(dark[: int(rb0), int(0.12 * W) : int(0.35 * W)].sum(axis=0)) + int(0.12 * W))

    # Leftward y-ticks: 4 upper ticks (7, 6, 5, 4) and 3 lower ticks (2, 1, 0)
    yt = [c[0] for c in clusters(np.where(dark[: int(rb1) + 5, v_col - 26 : v_col - 8].mean(axis=1) > 0.5)[0]) if (c[2] - c[1] + 1) <= 15]
    if len(yt) != 7:
        raise RuntimeError(f"Fig. 3B bar: expected 7 y-ticks (7,6,5,4,2,1,0), found {len(yt)}: {yt}")
    y_upper = verify_axis("linear", list(zip(yt[:4], [7.0, 6.0, 5.0, 4.0])), max_residual_px=1.5)
    y_lower = verify_axis("linear", list(zip(yt[4:], [2.0, 1.0, 0.0])), max_residual_px=1.5)
    r_break = 0.5 * (yt[3] + yt[4])

    # Three downward x-ticks for Wild-type (1), Src42A mutant (2), kkv mutant (3)
    xt = [c[0] for c in clusters(np.where(dark[rb1 + 8 : rb1 + 32, v_col + 20 :].mean(axis=0) > 0.5)[0] + v_col + 20) if (c[2] - c[1] + 1) <= 15]
    if len(xt) != 3:
        raise RuntimeError(f"Fig. 3B bar: expected 3 x-ticks, found {len(xt)}: {xt}")

    labels = ["Wild-type", "Src42A mutant", "kkv mutant"]
    out = []
    for idx, (xc, label) in enumerate(zip(xt, labels), start=1):
        xc_i = int(round(xc))
        # Bar top border measured outside error-bar cap width (xc-180..xc-130)
        top_rows = np.where((gray[int(yt[0]) - 20 : int(rb0) - 10, xc_i - 180 : xc_i - 130] < 120.0).mean(axis=1) > 0.7)[0] + (int(yt[0]) - 20)
        r_bar = float(top_rows.min()) + 0.5
        is_upper = r_bar < r_break
        y_ax = y_upper if is_upper else y_lower

        # Error-bar horizontal cap above r_bar (stay below axis break for lower bars, below top margin for upper bar)
        r_search_top = int(yt[0]) - 40 if is_upper else int(r_break) + 10
        cap_rows = np.where((gray[r_search_top : int(r_bar) - 10, xc_i - 45 : xc_i + 45] < 120.0).mean(axis=1) > 0.5)[0] + r_search_top
        if len(cap_rows) == 0:
            raise RuntimeError(f"Fig. 3B bar: could not locate error-bar cap for {label}")
        r_cap = float(cap_rows.mean()) + 0.5

        half_px = 2.0 + y_ax.residual_px
        y_val, lo, hi = value_interval(y_ax, r_bar, half_px)
        bar_hi = float(y_ax.to_value(r_cap))
        bar_hi_out = float(y_ax.to_value(r_cap - half_px))
        # Lower error-bar half hidden inside filled bar -> mirror of upper half
        bar_lo = 2.0 * y_val - bar_hi
        bar_lo_out = 2.0 * y_val - bar_hi_out

        out.append({
            "x": float(idx),
            "y": y_val,
            "lo": lo,
            "hi": hi,
            "bar_lo": bar_lo,
            "bar_hi": bar_hi,
            "bar_lo_outer": bar_lo_out,
            "bar_hi_outer": bar_hi_out,
            "px": [float(xc), r_bar],
            "px_lo": r_cap,
            "px_hi": r_bar,
            "note": f"{label} ({'upper' if is_upper else 'lower'} y-segment); lower error bar hidden inside bar (mirrored); cal_res={y_ax.residual_px:.2f}px",
        })
    return out
