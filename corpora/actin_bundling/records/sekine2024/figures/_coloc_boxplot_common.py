"""Shared reader for 7-protein Pearson's R colocalization boxplots (Fig. 3d Stage 15 and Supp. Fig. 5d Stage 16, 900 dpi)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage


def _split_runs(idx: np.ndarray, gap: int = 5) -> list[np.ndarray]:
    if len(idx) == 0:
        return []
    cuts = np.where(np.diff(idx) > gap)[0] + 1
    return [g for g in np.split(idx, cuts) if len(g) > 0]


def read_coloc_boxplots(
    image: np.ndarray,
    n_expected: int = 7,
    min_bar_w: int = 100,
    min_bar_h: int = 8,
    y_max_data: float = 1.02,
):
    """Calibrate y-axis (1.0 to -1.0) and read median bar + whisker/outlier range for all 7 protein boxplots."""
    img = image.astype(float)
    H, W, _ = img.shape
    g = img.mean(axis=2)

    hr = _split_runs(np.where((g < 120.0).sum(axis=1) > W * 0.45)[0], gap=10)[-1]
    vc = _split_runs(np.where((g < 120.0).sum(axis=0) > H * 0.40)[0], gap=10)[0]

    y_labels = [1.0, 0.5, 0.0, -0.5, -1.0]
    yp = (255.0 - g[:, max(0, int(vc.min()) - 15) : max(0, int(vc.min()) - 2)]).mean(axis=1)
    yp[int(hr.min()) + 5 :] = 0.0
    yt = [float(((gr + 0.5) * yp[gr]).sum() / yp[gr].sum()) for gr in _split_runs(np.where(yp > 28.0)[0], gap=4)]
    if len(yt) != len(y_labels):
        raise RuntimeError(f"Expected {len(y_labels)} y-ticks, found {len(yt)}: {yt}")

    by, ay = np.polyfit(yt, y_labels, 1)
    res_y = float(np.max(np.abs(np.array(yt) - (np.array(y_labels) - ay) / by)))
    if res_y > 2.0:
        raise RuntimeError(f"Y calibration residual too high: ry={res_y:.2f}px")

    r_top = max(10, int(round((y_max_data - ay) / by)))
    r_bot = min(int(hr.min()) - 5, int(round((-0.82 - ay) / by)))
    r_zero = (0.0 - ay) / by
    c_left = int(vc.max()) + 10
    c_right = W - 15

    # Detect the 7 x-axis tick centers below hr (exact stem column for whiskers & outlier dots)
    xp = (255.0 - g[int(hr.max()) + 2 : min(H, int(hr.max()) + 15), :]).mean(axis=0)
    xp[: int(vc.max()) + 20] = 0.0
    xt = [float(((gr + 0.5) * xp[gr]).sum() / xp[gr].sum()) for gr in _split_runs(np.where(xp > 25.0)[0], gap=4)]
    if len(xt) != n_expected:
        raise RuntimeError(f"Expected {n_expected} x-ticks, found {len(xt)}: {xt}")

    # Detect the thick horizontal median bars (opening kills thin box edges and jitter dots)
    sub = g[r_top:r_bot, c_left:c_right]
    thick_h = ndimage.binary_opening(sub < 165.0, structure=np.ones((min_bar_h, min_bar_w), dtype=bool))
    lab, n = ndimage.label(thick_h)

    boxes = []
    for k in range(1, n + 1):
        rr, cc = np.nonzero(lab == k)
        c_bar = 0.5 * (float(cc.min()) + float(cc.max())) + 0.5 + c_left
        r_med = float(rr.mean()) + 0.5 + r_top
        boxes.append((c_bar, r_med, int(cc.min()) + c_left, int(cc.max()) + c_left))
    boxes.sort(key=lambda b: b[0])
    if len(boxes) != n_expected:
        raise RuntimeError(f"Expected {n_expected} boxplot medians, found {len(boxes)}")

    out = []
    for i, (_c_bar, r_med, b_c0, b_c1) in enumerate(boxes):
        c_mid = xt[i]
        ci = int(round(c_mid))
        # Find box top/bottom borders at off-center column (avoiding central whisker stem)
        c_edge = max(c_left, ci - 60)
        edge_dark = g[r_top:r_bot, c_edge : c_edge + 15].mean(axis=1) < 160.0
        # Exclude the dashed y=0 reference line unless it touches a box edge
        edge_rows = [r + r_top for r in np.nonzero(edge_dark)[0] if abs((r + r_top) - r_zero) > 5.0]
        box_top_r = float(min(edge_rows)) if edge_rows else r_med
        box_bot_r = float(max(edge_rows)) if edge_rows else r_med

        # Find vertical whisker ends and any black outlier dots along the central stem ci ± 6
        col_dark = g[r_top:r_bot, max(0, ci - 6) : min(W, ci + 7)].min(axis=1) < 130.0
        runs = _split_runs(np.where(col_dark)[0], gap=4)
        # Keep vertical stem/box runs or outlier dots (height >= 7 px) and ignore lone y=0 horizontal dash
        valid = [gr for gr in runs if len(gr) >= 7 and not (len(gr) <= 8 and abs((gr.mean() + r_top) - r_zero) < 5.0)]
        r_hi_val = min(box_top_r, float(valid[0][0]) + 0.5 + r_top) if valid else box_top_r
        r_lo_val = max(box_bot_r, float(valid[-1][-1]) + 0.5 + r_top) if valid else box_bot_r
        # If the bottommost or topmost valid run is an isolated compact outlier dot (height <= 22px, separated from box), use its center!
        if valid and len(valid[-1]) <= 22 and (valid[-1][0] + r_top - box_bot_r) > 25:
            r_lo_val = float(valid[-1].mean()) + 0.5 + r_top
        if valid and len(valid[0]) <= 22 and (box_top_r - (valid[0][-1] + r_top)) > 25:
            r_hi_val = float(valid[0].mean()) + 0.5 + r_top

        y_med = ay + by * r_med
        bar_hi = max(y_med, ay + by * r_hi_val)
        bar_lo = min(y_med, ay + by * r_lo_val)
        err_px = 2.0 + res_y
        err_y = abs(by) * err_px
        bar_pad = abs(by) * (3.0 + res_y)
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
            "note": f"boxplot median + min/max range (group={i+1}, ry={res_y:.2f}px)",
        })
    return out
