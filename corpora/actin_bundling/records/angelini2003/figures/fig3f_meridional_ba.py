"""Calibrated reader for Fig. 3F (p. 3, main panel): Ba2+-condensed F-actin meridional SAXS cut I(q_z)/I_0."""

from __future__ import annotations

import numpy as np

from boxgraph.digitize import Axis


def _longest_run_len(arr: np.ndarray) -> int:
    if not arr.any():
        return 0
    padded = np.concatenate(([False], arr, [False]))
    diff = np.diff(padded.astype(int))
    starts = np.where(diff == 1)[0]
    ends = np.where(diff == -1)[0]
    return int((ends - starts).max())


def _cluster_1d(idxs: np.ndarray, gap: int = 4) -> list[float]:
    if len(idxs) == 0:
        return []
    groups = np.split(idxs, np.where(np.diff(idxs) > gap)[0] + 1)
    return [float(g.mean()) for g in groups if len(g) > 0]


def _calibrate_fig3f(gray: np.ndarray):
    H, W = gray.shape
    dark = gray < 170
    # Outer frame spans >70% of image width/height; inset bottom spans ~52% of width around middle-upper region
    h_outer = _cluster_1d(np.array([r for r in range(H) if _longest_run_len(dark[r]) > 0.70 * W]))
    v_outer = _cluster_1d(np.array([c for c in range(W) if _longest_run_len(dark[:, c]) > 0.70 * H]))
    r0, r1 = h_outer[0], h_outer[-1]
    c0, c1 = v_outer[0], v_outer[-1]

    # Locate inset bottom horizontal border (~r=882) so redness search stays strictly below the inset
    h_inset = _cluster_1d(np.array([r for r in range(int(r0) + 200, int(r1) - 400) if _longest_run_len(dark[r]) > 0.45 * W]))
    r_inset_bot = h_inset[-1] if h_inset else r0 + 0.42 * (r1 - r0)

    bt = gray[int(r1) + 10 : int(r1) + 34, max(0, int(c0) - 12) : min(W, int(c1) + 12)] < 170
    xt = [x + max(0, int(c0) - 12) for x in _cluster_1d(np.where(bt.mean(axis=0) > 0.30)[0])]
    x_vals = [0.05 + 0.025 * i for i in range(len(xt))]
    ax_x = Axis("linear", list(zip(xt, x_vals)))

    lt = gray[max(0, int(r0) - 10) : min(H, int(r1) + 10), max(0, int(c0) - 35) : max(1, int(c0) - 10)] < 175
    yt = [y + max(0, int(r0) - 10) for y in _cluster_1d(np.where(lt.mean(axis=1) > 0.25)[0])]
    major_y = max(
        yt,
        key=lambda y: int((gray[int(round(y)) - 2 : int(round(y)) + 3, max(0, int(c0) - 60) : int(c0) - 6].mean(axis=0) < 175).sum()),
    )
    above = [y for y in yt if y < major_y - 20]
    below = [y for y in yt if y > major_y + 20]
    y_pairs = [(major_y, 10.0)]
    if above:
        y_pairs.append((above[-1], 20.0))
    for k, y in enumerate(below):
        v = 9 - k
        if v >= 2:
            y_pairs.append((y, float(v)))
    ax_y = Axis("log10", y_pairs)
    return ax_x, ax_y, (c0, r0, c1, r1), r_inset_bot


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    arr = image.astype(float)
    gray = arr.mean(axis=2)
    redness = arr[:, :, 0] - 0.5 * (arr[:, :, 1] + arr[:, :, 2])
    ax_x, ax_y, (c0, r0, c1, r1), r_inset_bot = _calibrate_fig3f(gray)

    rs, re = int(r_inset_bot + 10), int(r1 - 6)
    px_quant = 4.0
    x_err_px = np.hypot(px_quant, ax_x.residual_px)
    y_cal_px = np.hypot(px_quant, ax_y.residual_px)

    def _col_row(ci: int):
        rr_list, ww_list = [], []
        for dc in range(-4, 5):
            cc = ci + dc
            if int(c0) + 4 <= cc <= int(c1) - 4:
                col = redness[rs:re, cc]
                m = float(col.max())
                if m > 12.0:
                    idx = np.where(col >= 0.45 * m)[0]
                    rr_list.append(idx + rs)
                    ww_list.append(col[idx])
        if not rr_list:
            return None, 0.0
        rr_all = np.concatenate(rr_list)
        ww_all = np.concatenate(ww_list)
        return float(np.sum(rr_all * ww_all) / np.sum(ww_all)), 0.5 * float(rr_all.max() - rr_all.min())

    out = []
    for q in at:
        c = ax_x.to_pixel(q)
        ci = int(round(c))
        row, half_thick = _col_row(ci)
        if row is None:
            out.append({"x": float(q), "y": None, "note": "Fig. 3F red curve not found at q_z"})
            continue
        rl, _ = _col_row(ci - 12)
        rr, _ = _col_row(ci + 12)
        slope_px = abs(rr - rl) / 24.0 if (rl is not None and rr is not None) else 0.0
        row_err = float(np.hypot(np.hypot(y_cal_px, 0.5 * half_thick), slope_px * x_err_px))
        y_val = ax_y.to_value(row)
        y_hi = ax_y.to_value(row - row_err)
        y_lo = ax_y.to_value(row + row_err)
        out.append({
            "x": float(q),
            "y": float(y_val),
            "lo": float(min(y_lo, y_hi)),
            "hi": float(max(y_lo, y_hi)),
            "px": [float(c), float(row)],
            "px_lo": float(row - row_err),
            "px_hi": float(row + row_err),
            "note": f"x_res={ax_x.residual_px:.2f}px, y_res={ax_y.residual_px:.2f}px",
        })
    return out
