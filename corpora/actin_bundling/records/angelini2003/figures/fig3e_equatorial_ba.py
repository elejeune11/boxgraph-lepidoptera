"""Calibrated reader for Fig. 3E (p. 3): Ba2+-condensed F-actin equatorial SAXS cut I(q_r)/I_0."""

from __future__ import annotations

import numpy as np
from scipy import ndimage

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


def _calibrate_fig3e(gray: np.ndarray):
    H, W = gray.shape
    dark = gray < 170
    h_rows = np.array([r for r in range(H) if _longest_run_len(dark[r]) > 0.70 * W])
    v_cols = np.array([c for c in range(W) if _longest_run_len(dark[:, c]) > 0.70 * H])
    hs = _cluster_1d(h_rows)
    vs = _cluster_1d(v_cols)
    r0, r1 = hs[0], hs[-1]
    c0, c1 = vs[0], vs[-1]

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
    return ax_x, ax_y, (c0, r0, c1, r1)


def _trace_grey_curve(gray: np.ndarray, box: tuple[float, float, float, float]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Erase pure-black arrows/letters (gray < 65, dilated) so only the medium-grey SAXS curve remains."""
    c0, r0, c1, r1 = box
    rs, re = int(r0 + 25), int(r1 - 10)
    cs, ce = int(c0 + 14), int(c1 - 14)
    sub = gray[rs:re, cs:ce].copy()
    annot = ndimage.binary_dilation(sub < 65.0, iterations=8)
    sub[annot] = 255.0

    cols, rows, half_w = [], [], []
    for j in range(sub.shape[1]):
        col = sub[:, j]
        dk = np.where(col < 228.0)[0]
        if len(dk) == 0:
            continue
        groups = np.split(dk, np.where(np.diff(dk) > 10)[0] + 1)
        valid = [gr for gr in groups if len(gr) >= 3 and col[gr].min() < 222.0]
        if not valid:
            continue
        gr = valid[-1]
        w = 255.0 - col[gr]
        rc = rs + float(np.sum(gr * w) / np.sum(w))
        cols.append(float(cs + j))
        rows.append(rc)
        half_w.append(0.5 * float(gr[-1] - gr[0]))
    return np.array(cols), np.array(rows), np.array(half_w)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    gray = image.astype(float).mean(axis=2)
    ax_x, ax_y, box = _calibrate_fig3e(gray)
    cols_arr, rows_arr, hw_arr = _trace_grey_curve(gray, box)

    px_quant = 4.0
    x_err_px = np.hypot(px_quant, ax_x.residual_px)
    y_cal_px = np.hypot(px_quant, ax_y.residual_px)

    out = []
    for q in at:
        c = ax_x.to_pixel(q)
        near = np.abs(cols_arr - c) <= 5.0
        if not near.any():
            out.append({"x": float(q), "y": None, "note": "Fig. 3E curve not found at q_r"})
            continue
        row = float(np.median(rows_arr[near]))
        half_thick = float(np.median(hw_arr[near]))
        left = (cols_arr >= c - 16) & (cols_arr <= c - 8)
        right = (cols_arr >= c + 8) & (cols_arr <= c + 16)
        slope_px = 0.0
        if left.any() and right.any():
            slope_px = abs(float(np.median(rows_arr[right]) - np.median(rows_arr[left]))) / 24.0
        elif left.any():
            slope_px = abs(row - float(np.median(rows_arr[left]))) / 12.0
        elif right.any():
            slope_px = abs(float(np.median(rows_arr[right])) - row) / 12.0

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
