"""Calibrated reader for Fig. 3C (p. 3): nonelectrostatically condensed F-actin meridional SAXS cut I(q_z)/I_0."""

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


def _calibrate_fig3c(gray: np.ndarray):
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
    x_vals = [0.025 + 0.025 * i for i in range(len(xt))]
    ax_x = Axis("linear", list(zip(xt, x_vals)))

    lt = gray[max(0, int(r0) - 10) : min(H, int(r1) + 10), max(0, int(c0) - 35) : max(1, int(c0) - 10)] < 175
    yt = [y + max(0, int(r0) - 10) for y in _cluster_1d(np.where(lt.mean(axis=1) > 0.25)[0])]
    majors = []
    for y in yt:
        ry = int(round(y))
        strip = gray[ry - 2 : ry + 3, max(0, int(c0) - 60) : int(c0) - 6].mean(axis=0) < 175
        if int(strip.sum()) > 24:
            majors.append(y)
    y_100, y_1 = majors[0], majors[-1]
    px_per_dec = 0.5 * (y_1 - y_100)
    cand_vals = (
        [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        + [1, 2, 3, 4, 5, 6, 7, 8, 9]
        + [10, 20, 30, 40, 50, 60, 70, 80, 90]
        + [100, 200, 300]
    )
    matched = []
    for y in yt:
        est = 10.0 ** ((y_1 - y) / px_per_dec)
        best_v = min(cand_vals, key=lambda v: abs(np.log10(v) - np.log10(est)))
        if abs(y - (y_1 - px_per_dec * np.log10(best_v))) < 10.0:
            matched.append((y, best_v))
    ax_y = Axis("log10", matched)
    return ax_x, ax_y, (c0, r0, c1, r1)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    arr = image.astype(float)
    gray = arr.mean(axis=2)
    redness = arr[:, :, 0] - 0.5 * (arr[:, :, 1] + arr[:, :, 2])
    ax_x, ax_y, (c0, r0, c1, r1) = _calibrate_fig3c(gray)

    rs, re = int(r0 + 10), int(r1 - 6)
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
            out.append({"x": float(q), "y": None, "note": "Fig. 3C red curve not found at q_z"})
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
