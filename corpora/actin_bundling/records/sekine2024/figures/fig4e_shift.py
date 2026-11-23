"""Digitize all 6 shift boxplots (1=1M_dx, 2=1M_dy, 3=2M_dx, 4=2M_dy, 5=4M_dx, 6=4M_dy) in Fig. 4e (PDF p. 7, 1200 dpi)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage


def _split_runs(idx: np.ndarray, gap: int = 5) -> list[np.ndarray]:
    if len(idx) == 0:
        return []
    cuts = np.where(np.diff(idx) > gap)[0] + 1
    return [g for g in np.split(idx, cuts) if len(g) > 0]


def read(image, at):
    img = image.astype(float)
    H, W, _ = img.shape
    g = img.mean(axis=2)

    # 3 vertical y-axes (one per subpanel: 1-M, 2-M, 4-M)
    v_axes = _split_runs(np.where((g < 120.0).sum(axis=0) > H * 0.50)[0], gap=10)
    if len(v_axes) != 3:
        raise RuntimeError(f"Expected 3 subpanel y-axes in Fig. 4e, found {len(v_axes)}")

    y_labels = [0.08, 0.06, 0.04, 0.02, 0.00]
    out = []
    cond_idx = 0

    for si, v in enumerate(v_axes):
        c_start = int(v.max()) + 10
        c_end = int(v_axes[si + 1].min()) - 180 if si + 1 < len(v_axes) else W - 20

        # Outward y-ticks directly left of this subpanel's y-axis
        yp = (255.0 - g[:, max(0, int(v.min()) - 18) : max(0, int(v.min()) - 2)]).mean(axis=1)
        yt = [float(((gr + 0.5) * yp[gr]).sum() / yp[gr].sum()) for gr in _split_runs(np.where(yp > 30.0)[0], gap=4)]
        if len(yt) != len(y_labels):
            raise RuntimeError(f"Panel {si}: expected 5 y-ticks, found {len(yt)}: {yt}")

        by, ay = np.polyfit(yt, y_labels, 1)
        res_y = float(np.max(np.abs(np.array(yt) - (np.array(y_labels) - ay) / by)))
        if res_y > 2.5:
            raise RuntimeError(f"Panel {si}: calibration residual too high: ry={res_y:.2f}px")

        # Find the horizontal x-axis line for this subpanel and its 2 outward x-ticks (dx, dy)
        hr_sub = _split_runs(np.where((g[:, c_start:c_end] < 120.0).sum(axis=1) > (c_end - c_start) * 0.50)[0], gap=10)[-1]
        xp = (255.0 - g[int(hr_sub.max()) + 2 : min(H, int(hr_sub.max()) + 18), c_start:c_end]).mean(axis=0)
        xt_local = [float(((gr + 0.5) * xp[gr]).sum() / xp[gr].sum()) for gr in _split_runs(np.where(xp > 30.0)[0], gap=4)]
        if len(xt_local) != 2:
            raise RuntimeError(f"Panel {si}: expected 2 x-ticks (dx, dy), found {len(xt_local)}")

        r_top_search = int(round((0.065 - ay) / by))
        r_bot_search = int(round((0.000 - ay) / by)) + 4
        r_box_top = int(round((0.038 - ay) / by))
        r_box_bot = int(round((0.003 - ay) / by))

        for bi, xl in enumerate(xt_local):
            cond_idx += 1
            cm = int(round(xl + c_start))
            # Inside each box (width ~350 px, centered on cm), find the thick dark median bar
            # Using the 25th percentile across 220 columns cleanly ignores pale translucent gray jitter dots on top of the bar!
            strip_g = g[r_box_top:r_box_bot, cm - 110 : cm + 111]
            row_p25 = np.percentile(strip_g, 25, axis=1)
            # Median bar is the thickest darkest horizontal band (< 125)
            med_runs = [gr for gr in _split_runs(np.where(row_p25 < 125.0)[0], gap=2) if len(gr) >= 5]
            if not med_runs:
                raise RuntimeError(f"Panel {si} box {bi}: median bar not found")
            best_med = max(med_runs, key=lambda gr: len(gr) * (255.0 - row_p25[gr].mean()))
            r_med = float(best_med.mean()) + 0.5 + r_box_top

            # Upper whisker above median: dark vertical stem (< 145 avoids pale translucent jitter cloud > 175)
            rm_i = int(round(r_med))
            w_up = g[r_top_search:rm_i, cm - 5 : cm + 6].min(axis=1) < 145.0
            up_runs = [gr for gr in _split_runs(np.where(w_up)[0], gap=5) if len(gr) >= 10]
            if not up_runs:
                raise RuntimeError(f"Panel {si} box {bi}: upper whisker not found")
            r_hi_val = float(up_runs[0][0]) + 0.5 + r_top_search

            # Lower whisker below median: dense translucent jitter dots near y=0 lighten the stem slightly (~150-160),
            # so use < 165 down to r_bot_search (strictly above bottom x-axis spine).
            w_dn = g[rm_i:r_bot_search, cm - 5 : cm + 6].min(axis=1) < 165.0
            dn_runs = [gr for gr in _split_runs(np.where(w_dn)[0], gap=8) if len(gr) >= 10]
            if not dn_runs:
                raise RuntimeError(f"Panel {si} box {bi}: lower whisker not found")
            r_lo_val = float(dn_runs[-1][-1]) + 0.5 + rm_i

            y_med = ay + by * r_med
            bar_hi = max(y_med, ay + by * r_hi_val)
            bar_lo = min(y_med, ay + by * r_lo_val)
            err_px = 2.0 + res_y
            err_y = abs(by) * err_px
            bar_pad = abs(by) * (3.0 + res_y)

            out.append({
                "x": float(cond_idx),
                "y": y_med,
                "lo": y_med - err_y,
                "hi": y_med + err_y,
                "bar_lo": bar_lo,
                "bar_hi": bar_hi,
                "bar_lo_outer": bar_lo - bar_pad,
                "bar_hi_outer": bar_hi + bar_pad,
                "px": [float(cm), r_med],
                "px_lo": r_hi_val,
                "px_hi": r_lo_val,
                "note": f"panel={si} box={bi} median + whisker range (ry={res_y:.2f}px)",
            })
    return out
