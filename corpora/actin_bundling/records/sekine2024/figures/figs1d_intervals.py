"""Digitize WT early stage 16 periodic cable interval boxplot (median + min/max range, n=294, N=8) in Supplementary Fig. 1d (PDF p. 15, 1200 dpi)."""

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

    vc = _split_runs(np.where((g < 150.0).sum(axis=0) > H * 0.45)[0], gap=10)[0]
    hr = _split_runs(np.where((g < 150.0).sum(axis=1) > W * 0.40)[0], gap=10)[-1]

    y_labels = [0.6, 0.5, 0.4, 0.3]
    yp = (255.0 - g[:, max(0, int(vc.min()) - 16) : max(0, int(vc.min()) - 2)]).mean(axis=1)
    yp[int(hr.min()) + 5 :] = 0.0
    yt = [float(((gr + 0.5) * yp[gr]).sum() / yp[gr].sum()) for gr in _split_runs(np.where(yp > 25.0)[0], gap=4)]
    if len(yt) != len(y_labels):
        raise RuntimeError(f"Expected 4 y-ticks in Supp. Fig. 1d, found {len(yt)}: {yt}")

    by, ay = np.polyfit(yt, y_labels, 1)
    res_y = float(np.max(np.abs(np.array(yt) - (np.array(y_labels) - ay) / by)))
    if res_y > 2.0:
        raise RuntimeError(f"Y calibration residual too high: ry={res_y:.2f}px")

    # Boxplot sits in interval [0.40, 0.50] um
    r_top = int(round((0.52 - ay) / by))
    r_bot = int(round((0.40 - ay) / by))
    c_left = int(vc.max()) + 30
    c_right = W - 30

    sub = g[r_top:r_bot, c_left:c_right]
    # Box top, thick median bar, and bottom edge each span ~400 px horizontally;
    # use dark pixel count per row so an overlapping gray dot does not split the median bar.
    row_cnt = (sub < 140.0).sum(axis=1)
    h_runs = _split_runs(np.where(row_cnt > 250)[0], gap=2)
    if len(h_runs) < 3:
        raise RuntimeError(f"Expected >=3 box horizontal lines in Supp. Fig. 1d, found {len(h_runs)}")
    # Median bar is the thickest dark horizontal run (6 px vs 4 px for box borders)
    med_run = max(h_runs, key=lambda gr: (len(gr), float((255.0 - sub[gr]).sum())))
    r_med = float(med_run.mean()) + 0.5 + r_top
    # Use top box border (no overlapping dot) to get exact horizontal span and center
    top_run = h_runs[0]
    box_cols = np.where(sub[top_run].min(axis=0) < 140.0)[0]
    c_box0 = int(box_cols.min()) + c_left
    c_box1 = int(box_cols.max()) + c_left
    c_mid = 0.5 * (float(c_box0) + float(c_box1)) + 0.5
    ci = int(round(c_mid))

    # Vertical whiskers along ci ± 6
    col_dark = g[r_top:r_bot, ci - 6 : ci + 7].min(axis=1) < 140.0
    w_rr = np.nonzero(col_dark)[0]
    r_wtop = float(w_rr.min()) + 0.5 + r_top
    r_wbot = float(w_rr.max()) + 0.5 + r_top

    # Also detect gray embryo dots inside the box horizontal span to include true plotted max/min (e.g. top dot ~0.476 um)
    dot_sub = g[r_top:r_bot, c_box0:c_box1]
    dots_m = ndimage.binary_opening((dot_sub > 145.0) & (dot_sub < 215.0), structure=np.ones((11, 11), dtype=bool))
    dlab, dn = ndimage.label(dots_m)
    dot_rs = [float(np.nonzero(dlab == k)[0].mean()) + 0.5 + r_top for k in range(1, dn + 1) if (dlab == k).sum() >= 150]
    r_hi_val = min(r_wtop, min(dot_rs)) if dot_rs else r_wtop
    r_lo_val = max(r_wbot, max(dot_rs)) if dot_rs else r_wbot

    y_med = ay + by * r_med
    bar_hi = ay + by * r_hi_val
    bar_lo = ay + by * r_lo_val
    err_px = 2.0 + res_y
    err_y = abs(by) * err_px
    bar_pad = abs(by) * (3.0 + res_y)

    return [{
        "x": 1.0,
        "y": y_med,
        "lo": y_med - err_y,
        "hi": y_med + err_y,
        "bar_lo": bar_lo,
        "bar_hi": bar_hi,
        "bar_lo_outer": bar_lo - bar_pad,
        "bar_hi_outer": bar_hi + bar_pad,
        "px": [c_mid, r_med],
        "px_lo": r_hi_val,
        "px_hi": r_lo_val,
        "note": f"WT cable interval median + min/max range (dots={len(dot_rs)}, ry={res_y:.2f}px)",
    }]
