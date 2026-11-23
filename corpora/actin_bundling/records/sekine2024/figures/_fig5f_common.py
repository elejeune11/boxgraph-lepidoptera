"""Shared reader for Fig. 5f (Early stage 16 PIV mean displacement boxplots across 4 genotypes, PDF p. 9, 1200 dpi)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage


def _split_runs(idx: np.ndarray, gap: int = 5) -> list[np.ndarray]:
    if len(idx) == 0:
        return []
    cuts = np.where(np.diff(idx) > gap)[0] + 1
    return [g for g in np.split(idx, cuts) if len(g) > 0]


def read_fig5f(image: np.ndarray, component: str):
    """Read Vx (coral-pink boxes) or Vy (cyan boxes) medians + min-max range across 1=Control, 2=DAAM RNAi, 3=Zasp52 RNAi, 4=myosin II DN."""
    img = image.astype(float)
    H, W, _ = img.shape
    g = img.mean(axis=2)

    hr = _split_runs(np.where((g < 150.0).sum(axis=1) > W * 0.45)[0], gap=10)[-1]
    vc = _split_runs(np.where((g < 150.0).sum(axis=0) > H * 0.45)[0], gap=10)[0]

    y_labels = [40.0, 30.0, 20.0, 10.0]
    yp = (255.0 - g[:, max(0, int(vc.min()) - 25) : max(0, int(vc.min()) - 4)]).mean(axis=1)
    yp[int(hr.min()) + 5 :] = 0.0
    yt = [float(((gr + 0.5) * yp[gr]).sum() / yp[gr].sum()) for gr in _split_runs(np.where(yp > 30.0)[0], gap=4)]
    if len(yt) != len(y_labels):
        raise RuntimeError(f"Expected 4 y-ticks in Fig. 5f, found {len(yt)}: {yt}")

    by, ay = np.polyfit(yt, y_labels, 1)
    res_y = float(np.max(np.abs(np.array(yt) - (np.array(y_labels) - ay) / by)))
    if res_y > 2.5:
        raise RuntimeError(f"Y calibration residual too high: ry={res_y:.2f}px")

    r_top_lim = int(round((39.1 - ay) / by))
    r_bot_lim = int(round((9.5 - ay) / by))
    c_left = int(vc.max()) + 10
    c_right = W - 150  # exclude Vx/Vy legend on far right

    sub = img[r_top_lim:r_bot_lim, c_left:c_right]
    R, G, B = sub[:, :, 0], sub[:, :, 1], sub[:, :, 2]

    if component == "vx":
        bmask = (R - G > 35.0) & (R - B > 35.0) & (R > 175.0)
    elif component == "vy":
        bmask = (G - R > 45.0) & (B - R > 45.0) & (B > 155.0)
    else:
        raise ValueError(f"Unknown component: {component}")

    # Open by (12, 70) to remove scatter dots (diam ~30px) while keeping colored box fills, then close vertically across median bar
    boxes_only = ndimage.binary_opening(bmask, structure=np.ones((12, 70), dtype=bool))
    joined = ndimage.binary_closing(boxes_only, structure=np.ones((40, 1), dtype=bool))
    lab, n = ndimage.label(joined)

    boxes = []
    for k in range(1, n + 1):
        rr, cc = np.nonzero(lab == k)
        if len(rr) < 1200:
            continue
        c0b = int(cc.min()) + c_left
        c1b = int(cc.max()) + c_left
        r0b = int(rr.min()) + r_top_lim
        r1b = int(rr.max()) + r_top_lim
        boxes.append((c0b, c1b, r0b, r1b))
    boxes.sort(key=lambda b: b[0])
    if len(boxes) != 4:
        raise RuntimeError(f"Expected 4 {component} boxes in Fig. 5f, found {len(boxes)}")

    out = []
    for i, (c0b, c1b, r0b, r1b) in enumerate(boxes):
        cm = 0.5 * (c0b + c1b) + 0.5
        ci = int(round(cm))

        # Thick dark median bar inside box (use 20th percentile across box width to ignore central colored dots!)
        r_s0 = max(r_top_lim, r0b - 18)
        r_s1 = min(r_bot_lim, r1b + 18)
        row_p20 = np.percentile(g[r_s0:r_s1, c0b + 12 : c1b - 12], 20, axis=1)
        # Median bar is the thickest darkest band (< 105, height >= 7)
        med_runs = [gr for gr in _split_runs(np.where(row_p20 < 105.0)[0], gap=2) if len(gr) >= 7]
        if not med_runs:
            raise RuntimeError(f"Fig. 5f {component} group {i+1}: median bar not found")
        best_med = max(med_runs, key=lambda gr: len(gr) * (255.0 - row_p20[gr].mean()))
        r_med = float(best_med.mean()) + 0.5 + r_s0

        # Full range (min-max plotted data dots of same color along central column ci ± 10; ignores black significance bars!)
        strip = img[r_top_lim:r_bot_lim, max(0, ci - 10) : min(W, ci + 11)]
        Rs, Gs, Bs = strip[:, :, 0], strip[:, :, 1], strip[:, :, 2]
        if component == "vx":
            col_color = (Rs - Gs > 45.0) & (Rs - Bs > 45.0) & (Rs > 175.0)
        else:
            col_color = (Gs - Rs > 45.0) & (Bs - Rs > 45.0) & (Bs > 140.0)
        d_runs = _split_runs(np.where(col_color.any(axis=1))[0], gap=4)
        valid_runs = [gr for gr in d_runs if len(gr) >= 12]
        if not valid_runs:
            raise RuntimeError(f"Fig. 5f {component} group {i+1}: range not found")

        # Extreme points are circular dots of radius ~13 px (or box edge if a dot touches box)
        top_run = valid_runs[0]
        bot_run = valid_runs[-1]
        r_hi_val = (float(top_run.mean()) + 0.5 if len(top_run) <= 36 else float(top_run[0]) + 13.0) + r_top_lim
        r_lo_val = (float(bot_run.mean()) + 0.5 if len(bot_run) <= 36 else float(bot_run[-1]) - 13.0) + r_top_lim

        y_med = ay + by * r_med
        bar_hi = max(y_med, ay + by * r_hi_val)
        bar_lo = min(y_med, ay + by * r_lo_val)
        err_px = 2.0 + res_y
        err_y = abs(by) * err_px
        bar_pad = abs(by) * (4.0 + res_y)

        out.append({
            "x": float(i + 1),
            "y": y_med,
            "lo": y_med - err_y,
            "hi": y_med + err_y,
            "bar_lo": bar_lo,
            "bar_hi": bar_hi,
            "bar_lo_outer": bar_lo - bar_pad,
            "bar_hi_outer": bar_hi + bar_pad,
            "px": [cm, r_med],
            "px_lo": r_hi_val,
            "px_hi": r_lo_val,
            "note": f"PIV {component} boxplot median + min/max range (group={i+1}, ry={res_y:.2f}px)",
        })
    return out
