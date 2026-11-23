"""Shared reader for black-diamond mean fusion durations in grouped boxplots (Fig. 2b, Fig. 3e, 900 dpi)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage


def _split_runs(idx: np.ndarray, gap: int = 5) -> list[np.ndarray]:
    if len(idx) == 0:
        return []
    cuts = np.where(np.diff(idx) > gap)[0] + 1
    return [g for g in np.split(idx, cuts) if len(g) > 0]


def read_fusion_means(
    image: np.ndarray,
    group_x_values: list[float],
    angle_slot: int,
    y_labels_top_to_bottom: list[float] = (9.0, 6.0, 3.0, 0.0),
):
    """Detect black diamond mean markers for one angle bin (0=60-90 deg, 1=30-60 deg, 2=0-30 deg) across groups."""
    img = image.astype(float)
    H, W, _ = img.shape
    gray = img.mean(axis=2)

    hr = _split_runs(np.where((gray < 100.0).sum(axis=1) > W * 0.50)[0], gap=10)[-1]
    vc = _split_runs(np.where((gray < 100.0).sum(axis=0) > H * 0.50)[0], gap=10)[0]

    # Outward y-ticks
    yp = (255.0 - gray[:, max(0, int(vc.min()) - 18) : max(0, int(vc.min()) - 3)]).mean(axis=1)
    yp[int(hr.min()) + 5 :] = 0.0
    yt = [float(((g + 0.5) * yp[g]).sum() / yp[g].sum()) for g in _split_runs(np.where(yp > 50.0)[0], gap=4)]
    if len(yt) != len(y_labels_top_to_bottom):
        raise RuntimeError(f"Expected {len(y_labels_top_to_bottom)} y-ticks, found {len(yt)}: {yt}")

    by, ay = np.polyfit(yt, list(y_labels_top_to_bottom), 1)
    res_y = float(np.max(np.abs(np.array(yt) - (np.array(y_labels_top_to_bottom) - ay) / by)))
    if res_y > 2.0:
        raise RuntimeError(f"Y calibration residual too high: ry={res_y:.2f}px")

    # Pure black diamonds strictly inside duration band [0.9, 6.0] frames (above n=... labels, below p-values)
    r_top = max(10, int(round((6.0 - ay) / by)))
    r_bot = min(int(hr.min()) - 5, int(round((0.9 - ay) / by)))
    c_left = int(vc.max()) + 15
    c_right = W - 60

    sub = gray[r_top:r_bot, c_left:c_right]
    # Pure black (< 25) opened by 15x15 square removes all dark-gray box borders and whiskers, leaving only solid black diamonds
    blk = ndimage.binary_opening(sub < 25.0, structure=np.ones((15, 15), dtype=bool))
    lab, n = ndimage.label(blk)

    diams = []
    for k in range(1, n + 1):
        rr, cc = np.nonzero(lab == k)
        if len(rr) < 200:
            continue
        c_abs = float(cc.mean()) + 0.5 + c_left
        r_abs = float(rr.mean()) + 0.5 + r_top
        diams.append((c_abs, r_abs, len(rr)))
    diams.sort(key=lambda d: d[0])

    n_expected = 3 * len(group_x_values)
    if len(diams) != n_expected:
        raise RuntimeError(f"Expected {n_expected} black diamonds (3 per group), found {len(diams)}")

    out = []
    for gi, xv in enumerate(group_x_values):
        c_abs, r_abs, area = diams[3 * gi + angle_slot]
        y_val = ay + by * r_abs
        err_px = 2.0 + res_y
        err_y = abs(by) * err_px
        out.append({
            "x": float(xv),
            "y": y_val,
            "lo": y_val - err_y,
            "hi": y_val + err_y,
            "px": [c_abs, r_abs],
            "px_lo": r_abs - err_px,
            "px_hi": r_abs + err_px,
            "note": f"black diamond mean (slot={angle_slot}, area={area}px, ry={res_y:.2f}px)",
        })
    return out
