"""Shared reader for SI Fig. 5A-C adult H. sara beeswarm + boxplots (PDF p. 28, 1200 dpi)."""

from __future__ import annotations

import numpy as np
from boxgraph import digitize as D
from _common import find_runs, run_centers, y_interval


def read_figs5_boxplots(image: np.ndarray, y_labels: list[float],
                        inset_row_max: int = 1100, inset_col_min: int = 4480) -> list[dict]:
    """Calibrate L-axes and extract median, IQR (Q1, Q3), whiskers, and full min-max data range
    for all 4 groups: 1=black cover, 2=black ground, 3=blue cover, 4=blue ground."""
    gray = image.astype(float).mean(axis=2)
    dark = gray < 100.0
    H, W = dark.shape

    c_frac = dark.mean(axis=0)
    c_frac[:200] = 0.0
    r_frac = dark.mean(axis=1)
    r_frac[:200] = 0.0
    c0 = run_centers(np.where(c_frac > 0.50)[0], gap=5)[0]
    r1 = run_centers(np.where(r_frac > 0.50)[0], gap=5)[0]

    # Ignore top panel border (row < 100) when detecting y-ticks
    yt = run_centers(np.where(dark[100:int(r1), int(round(c0)) - 25])[0] + 100, gap=5)
    xt = run_centers(np.where(dark[int(round(r1)) + 25, int(round(c0)):])[0] + int(round(c0)), gap=5)
    if len(yt) != len(y_labels) or len(xt) != 2:
        raise RuntimeError(f"SI Fig. 5 expected {len(y_labels)} y-ticks and 2 x-ticks, got yt={yt}, xt={xt}")
    ya = D.Axis("linear", list(zip(yt, y_labels)))

    # Box borders and whiskers are dark ink strictly inside plot (below top border row 60, left of inset)
    box_dark = gray < 85.0
    box_dark[:60, :] = False
    box_dark[int(r1) - 25:, :] = False
    box_dark[:, :int(c0) + 50] = False
    box_dark[:inset_row_max, inset_col_min:] = False

    # Colored data dots (cyan cover, gold ground) for full min-max range including outliers
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    colored = ((B - R > 20.0) & (G - R > 20.0)) | ((R - B > 30.0) & (G - B > 20.0))
    colored[:30, :] = False
    colored[int(r1) - 8:, :] = False
    colored[:, :int(c0) + 40] = False
    colored[:inset_row_max, inset_col_min:] = False

    dx = xt[1] - xt[0]
    centers_est = [xt[0] - 0.25 * dx, xt[0] + 0.25 * dx, xt[1] - 0.25 * dx, xt[1] + 0.25 * dx]
    names = ["black cover", "black ground", "blue cover", "blue ground"]

    out = []
    for gi, (ce, gname) in enumerate(zip(centers_est, names), start=1):
        cl = int(round(ce - 0.18 * dx))
        ch = int(round(ce + 0.18 * dx))
        sub_d = box_dark[:int(r1) - 25, cl:ch]
        row_counts = sub_d.sum(axis=1)
        h_runs = find_runs(np.where(row_counts >= 90)[0], gap=4)
        if len(h_runs) != 3:
            raise RuntimeError(f"Group {gi} ({gname}): expected 3 box horizontal edges (Q3, med, Q1), got {h_runs}")
        r_q3, r_med, r_q1 = h_runs[0][2], h_runs[1][2], h_runs[2][2]
        med_thick = 0.5 * (h_runs[1][1] - h_runs[1][0] + 1)

        # Exact horizontal box center
        cols_box = np.where(sub_d[int(round(r_med))])[0] + cl
        xc = float(cols_box.mean() + 0.5)

        # Whisker vertical stem along xc +- 5 px
        stem_rows = np.where(box_dark[:int(r1) - 25, int(round(xc)) - 5:int(round(xc)) + 6].any(axis=1))[0]
        r_w_top = float(stem_rows.min() + 0.5)
        r_w_bot = float(stem_rows.max() + 0.5)

        # Full data min-max from colored dots in this group's beeswarm column
        dl = max(int(c0) + 40, int(round(xc - 0.22 * dx)))
        dh = min(W - 5, int(round(xc + 0.22 * dx)))
        # Keep strictly left of inset when above inset_row_max
        dot_rows = np.where(colored[:int(r1) - 8, dl:dh].any(axis=1))[0]
        # Top/bottom of extreme dot circles have ~15 px radius; adjust by half a dot height (~14 px) to dot center
        r_d_top = float(dot_rows.min() + 14.0)
        r_d_bot = float(dot_rows.max() - 14.0)

        r_ext_top = min(r_w_top, r_d_top)
        r_ext_bot = max(r_w_bot, r_d_bot)

        y_med, lo, hi, _, _ = y_interval(ya, r_med, max(med_thick, 2.0))
        q3 = float(ya.to_value(r_q3))
        q1 = float(ya.to_value(r_q1))
        w_hi = float(ya.to_value(r_w_top))
        w_lo = float(ya.to_value(r_w_bot))
        d_hi = float(ya.to_value(r_ext_top))
        d_lo = float(ya.to_value(r_ext_bot))

        h_px = 4.0 + max(ya.residual_px, 0.5)
        bar_hi_outer = float(ya.to_value(r_ext_top - h_px))
        bar_lo_outer = float(ya.to_value(r_ext_bot + h_px))

        out.append({
            "x": float(gi),
            "y": float(y_med), "lo": float(lo), "hi": float(hi),
            "bar_lo": float(min(d_lo, d_hi)), "bar_hi": float(max(d_lo, d_hi)),
            "bar_lo_outer": float(min(bar_lo_outer, bar_hi_outer)),
            "bar_hi_outer": float(max(bar_lo_outer, bar_hi_outer)),
            "px": [xc, float(r_med)],
            "px_lo": float(r_ext_top), "px_hi": float(r_ext_bot),
            "note": (f"{gname}: med={y_med:.4g}, IQR=[{q1:.4g},{q3:.4g}], "
                     f"whisk=[{w_lo:.4g},{w_hi:.4g}], range=[{d_lo:.4g},{d_hi:.4g}], ry={ya.residual_px:.2f} px"),
        })
    return out
