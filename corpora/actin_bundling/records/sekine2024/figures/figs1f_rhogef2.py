"""Digitize Junctional RhoGEF2 intensity ratio boxplots (1=Stage 15, 2=Stage 16, n=5 each) in Supplementary Fig. 1f (PDF p. 15, 1200 dpi)."""

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

    y_labels = [1.5, 1.0, 0.5, 0.0]
    yp = (255.0 - g[:, max(0, int(vc.min()) - 16) : max(0, int(vc.min()) - 2)]).mean(axis=1)
    yp[int(hr.min()) + 5 :] = 0.0
    yt = [float(((gr + 0.5) * yp[gr]).sum() / yp[gr].sum()) for gr in _split_runs(np.where(yp > 25.0)[0], gap=4)]
    if len(yt) != len(y_labels):
        raise RuntimeError(f"Expected 4 y-ticks in Supp. Fig. 1f, found {len(yt)}: {yt}")

    by, ay = np.polyfit(yt, y_labels, 1)
    res_y = float(np.max(np.abs(np.array(yt) - (np.array(y_labels) - ay) / by)))
    if res_y > 2.0:
        raise RuntimeError(f"Y calibration residual too high: ry={res_y:.2f}px")

    # Detect the 2 bottom x-axis ticks (Stage 15, Stage 16) strictly right of the y-axis
    xp = (255.0 - g[int(hr.max()) + 2 : min(H, int(hr.max()) + 15), :]).mean(axis=0)
    xp[: int(vc.max()) + 20] = 0.0
    xt = [float(((gr + 0.5) * xp[gr]).sum() / xp[gr].sum()) for gr in _split_runs(np.where(xp > 25.0)[0], gap=4)]
    if len(xt) != 2:
        raise RuntimeError(f"Expected 2 x-ticks (Stage 15, Stage 16) in Supp. Fig. 1f, found {len(xt)}: {xt}")

    r_top = int(round((1.45 - ay) / by))
    r_bot = int(round((0.15 - ay) / by))

    out = []
    for i, c_mid in enumerate(xt):
        ci = int(round(c_mid))
        c0 = max(int(vc.max()) + 5, ci - 200)
        c1 = min(W - 5, ci + 200)
        sub = g[r_top:r_bot, c0:c1]
        row_cnt = (sub < 140.0).sum(axis=1)
        h_runs = _split_runs(np.where(row_cnt > 180)[0], gap=2)
        if len(h_runs) < 2:
            raise RuntimeError(f"Expected >=2 horizontal box lines for group {i+1} in Supp. Fig. 1f, found {len(h_runs)}")
        med_run = max(h_runs, key=lambda gr: (len(gr), float((255.0 - sub[gr]).sum())))
        r_med = float(med_run.mean()) + 0.5 + r_top

        # Find extreme top/bottom plotted data points (gray dots + vertical whiskers) along ci ± 10
        strip = g[r_top:r_bot, ci - 10 : ci + 11]
        # Topmost and bottommost gray dot centers (radius ~9 px) or whisker ends
        col_any = strip.min(axis=1) < 195.0
        w_rr = np.nonzero(col_any)[0]
        # Extreme top dot sits on top of upper whisker; dot radius is ~8px so its center is ~8px below topmost gray pixel
        # Verify by checking gray dot mask (125 < g < 200) near top and bottom
        top_patch = strip[max(0, w_rr.min() - 2) : min(strip.shape[0], w_rr.min() + 20), :]
        bot_patch = strip[max(0, w_rr.max() - 19) : min(strip.shape[0], w_rr.max() + 3), :]
        top_gray = np.nonzero((top_patch > 120.0) & (top_patch < 200.0))[0]
        bot_gray = np.nonzero((bot_patch > 120.0) & (bot_patch < 200.0))[0]
        if len(top_gray) >= 30:
            r_hi_val = float(top_gray.mean()) + max(0, w_rr.min() - 2) + 0.5 + r_top
        else:
            r_hi_val = float(w_rr.min()) + 0.5 + r_top
        if len(bot_gray) >= 30:
            r_lo_val = float(bot_gray.mean()) + max(0, w_rr.max() - 19) + 0.5 + r_top
        else:
            r_lo_val = float(w_rr.max()) + 0.5 + r_top

        y_med = ay + by * r_med
        bar_hi = max(y_med, ay + by * r_hi_val)
        bar_lo = min(y_med, ay + by * r_lo_val)
        err_px = 2.0 + res_y
        err_y = abs(by) * err_px
        bar_pad = abs(by) * (3.5 + res_y)
        out.append({
            "x": float(i + 1),
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
            "note": f"RhoGEF2 intensity ratio median + min/max range (stage={15+i}, ry={res_y:.2f}px)",
        })
    return out
