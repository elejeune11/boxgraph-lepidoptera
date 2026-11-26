"""Calibrated reader for Fig. 3C (p. 3): resolved meridional 13/6 layer-line peaks (d, f) of nonelectrostatically condensed F-actin."""

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


def _calibrate_x(gray: np.ndarray):
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
    return Axis("linear", list(zip(xt, x_vals))), (c0, r0, c1, r1)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    arr = image.astype(float)
    gray = arr.mean(axis=2)
    redness = arr[:, :, 0] - 0.5 * (arr[:, :, 1] + arr[:, :, 2])
    ax_x, (c0, r0, c1, r1) = _calibrate_x(gray)

    rs, re = int(r0 + 10), int(r1 - 6)
    cs, ce = int(c0 + 6), int(c1 - 6)
    cols, rows = [], []
    for c in range(cs, ce):
        col = redness[rs:re, c]
        m = float(col.max())
        if m > 12.0:
            w = np.maximum(col - 0.5 * m, 0.0)
            cols.append(float(c))
            rows.append(rs + float(np.sum(np.arange(len(col)) * w) / np.sum(w)))

    cols_arr = np.array(cols)
    rows_arr = np.array(rows)
    qs_arr = np.array([ax_x.to_value(c) for c in cols_arr])

    # Resolved local maxima: peak (d) 6th layer line (~0.111) and peak (f) 13th layer line (~0.225).
    # Peak (e) 7th layer line (~0.123) appears as a flat right shoulder on peak (d) and is captured in fig3c_meridional_nonelec.
    windows = {
        1: (0.104, 0.118, "meridional 6th layer-line peak (d)"),
        2: (0.215, 0.238, "meridional 13th layer-line peak (f)"),
    }
    cal_err_px = np.hypot(4.0, ax_x.residual_px)
    dq_per_px = abs(ax_x.b)

    out = []
    for idx_val in (at or [1.0, 2.0]):
        idx = int(round(idx_val))
        qlo, qhi, desc = windows[idx]
        m = (qs_arr >= qlo) & (qs_arr <= qhi)
        sub_c, sub_r = cols_arr[m], rows_arr[m]
        r_min = float(sub_r.min())
        crest = sub_r <= r_min + 8.0
        c_peak = float(sub_c[crest].mean())
        r_peak = float(sub_r[crest].mean())
        q_peak = float(ax_x.to_value(c_peak))
        crest_half_px = 0.5 * float(sub_c[crest].max() - sub_c[crest].min())
        c_err_px = float(np.hypot(cal_err_px, max(4.0, 0.5 * crest_half_px)))
        dq = c_err_px * dq_per_px
        out.append({
            "x": float(idx),
            "y": q_peak,
            "lo": float(q_peak - dq),
            "hi": float(q_peak + dq),
            "px": [c_peak, r_peak],
            "note": f"{desc}; x_res={ax_x.residual_px:.2f}px, crest_w={2*crest_half_px:.0f}px",
        })
    return out
