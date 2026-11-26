"""Calibrated reader for Fig. 3F Inset (p. 3): resolved meridional peaks i (CDW), j (l=17), k (l=19), and l (mosaic h)."""

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


def _calibrate_inset_x(gray: np.ndarray):
    H, W = gray.shape
    dark = gray < 180
    h_rows = _cluster_1d(np.array([r for r in range(H) if _longest_run_len(dark[r]) > 0.70 * W]))
    r0, r1 = h_rows[0], h_rows[-1]
    # Right border is solid dark; left border is slightly lighter grey
    v_right = _cluster_1d(np.array([c for c in range(W // 2, W) if _longest_run_len(dark[:, c]) > 0.65 * H]))
    c1 = v_right[-1]
    col_means = gray[int(r0) + 20 : int(r1) - 20, :150].mean(axis=0)
    c0 = float(np.argmin(col_means))

    # Bottom outward ticks just below r1: 5 evenly spaced ticks 0.10, 0.11, 0.12, 0.13, 0.14 A^-1
    prof = gray[int(r1) + 6 : int(r1) + 20, :].mean(axis=0)
    dk = np.where(prof < 205)[0]
    xt = [x for x in _cluster_1d(dk) if c0 + 5 <= x <= c1 - 15]
    x_vals = [0.10 + 0.01 * i for i in range(len(xt))]
    return Axis("linear", list(zip(xt, x_vals))), (c0, r0, c1, r1)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    arr = image.astype(float)
    gray = arr.mean(axis=2)
    redness = arr[:, :, 0] - 0.5 * (arr[:, :, 1] + arr[:, :, 2])
    ax_x, (c0, r0, c1, r1) = _calibrate_inset_x(gray)

    rs, re = int(r0 + 12), int(r1 - 8)
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

    # Four resolved peaks in Fig. 3F Inset: i (CDW), j (l=17), k (l=19), l (mosaic h)
    windows = {
        1: (0.1025, 0.1085, "sharp 1D counterion CDW peak (i)"),
        2: (0.1110, 0.1175, "36/17 l=17 layer-line peak (j)"),
        3: (0.1185, 0.1235, "36/17 l=19 layer-line peak (k)"),
        4: (0.1320, 0.1390, "mosaic-smeared equatorial peak (l)"),
    }
    cal_err_px = np.hypot(4.0, ax_x.residual_px)
    dq_per_px = abs(ax_x.b)

    out = []
    for idx_val in (at or [1.0, 2.0, 3.0, 4.0]):
        idx = int(round(idx_val))
        qlo, qhi, desc = windows[idx]
        m = (qs_arr >= qlo) & (qs_arr <= qhi)
        sub_c, sub_r = cols_arr[m], rows_arr[m]
        r_min = float(sub_r.min())
        crest = sub_r <= r_min + 10.0
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
