"""Shared calibrated stacked-panel time-series reader for CytoSim actin-covered surface area fraction
(Figs. 4J, 5C, S6B, S6F, S7D).

Isolates solid mean curves by strict RGB distance + saturation/luminance filtering so semi-transparent
SD shading bands are completely rejected.
"""

from __future__ import annotations

import numpy as np
import _common as C
from _bars_stacked import find_stacked_boxes


def read_ts_stacked(
    img: np.ndarray,
    n_rows: int,
    row_idx: int,
    color_rgb: tuple[float, float, float],
    color_tol: float,
    at: list[float],
    x_tick_vals: list[float] = (0.0, 200.0, 400.0, 600.0),
    y_tick_vals: list[float] = (0.6, 0.3, 0.0),
    max_lum: float | None = None,
    min_sat: float | None = None,
    max_sat: float | None = None,
    max_res_px: float = 4.0,
    search_half_cols: int = 3,
    label: str = "",
) -> list[dict]:
    """Read one solid mean time-series curve at timestamps `at` from row `row_idx` of `n_rows` stacked boxes."""
    boxes = find_stacked_boxes(img, n_rows)
    b = boxes[row_idx]
    c0, r0, c1, r1 = b["box"]

    yt = C.match_regular_ticks(b["y_ticks"], len(y_tick_vals), tol=6.0)
    ya = C.make_axis("linear", list(zip(yt, y_tick_vals)), max_res_px=max_res_px)

    xt_raw = boxes[-1]["x_ticks"]
    xt = C.match_regular_ticks(xt_raw, len(x_tick_vals), tol=6.0)
    xa = C.make_axis("linear", list(zip(xt, x_tick_vals)), max_res_px=max_res_px)

    sub = img[r0 + 2 : r1 - 1, c0 + 2 : c1 - 1].astype(np.float64)
    target = np.array(color_rgb, dtype=np.float64)
    dist = np.sqrt(((sub - target[None, None, :]) ** 2).sum(axis=2))
    mask = dist <= color_tol
    if max_lum is not None:
        mask &= sub.mean(axis=2) <= max_lum
    if min_sat is not None:
        sat = sub.max(axis=2) - sub.min(axis=2)
        mask &= sat >= min_sat
    if max_sat is not None:
        sat = sub.max(axis=2) - sub.min(axis=2)
        mask &= sat <= max_sat

    H_s, W_s = mask.shape
    out = []
    for xv in at:
        pc = xa.to_pixel(float(xv))
        c_rel = int(round(pc)) - (c0 + 2)
        if not (0 <= c_rel < W_s):
            out.append({"x": float(xv), "y": None, "note": f"{label}: x={xv} outside plot box"})
            continue

        # Gather matched curve pixels in narrow column window around c_rel
        cl = max(0, c_rel - search_half_cols)
        cr = min(W_s - 1, c_rel + search_half_cols)
        # Prefer exact column +/- 1 first
        ys_tight = np.nonzero(mask[:, max(0, c_rel - 1) : min(W_s, c_rel + 2)].any(axis=1))[0]
        ys = ys_tight if ys_tight.size > 0 else np.nonzero(mask[:, cl : cr + 1].any(axis=1))[0]
        if ys.size == 0:
            # Slightly widen search to +/- 6 cols in case another curve briefly crosses at exact x
            ys = np.nonzero(mask[:, max(0, c_rel - 6) : min(W_s, c_rel + 7)].any(axis=1))[0]
        if ys.size == 0:
            out.append({"x": float(xv), "y": None, "note": f"{label}: no solid pixels of {color_rgb} near t={xv}"})
            continue

        # Cluster vertical runs; pick the thickest contiguous run (solid line stroke, rejecting stray 1-px noise)
        cls = C.clusters_1d(ys, gap=3)
        best = max(cls, key=lambda t: t[2] - t[1] + 1)
        row_abs = best[0] + (r0 + 2)
        half_px = max(1.5, 0.5 * (best[2] - best[1] + 1))

        # Local slope estimate across +/- 6 columns
        slope = 0.0
        yl = np.nonzero(mask[:, max(0, c_rel - 6)])[0]
        yr = np.nonzero(mask[:, min(W_s - 1, c_rel + 6)])[0]
        if yl.size > 0 and yr.size > 0:
            c_l = min(C.clusters_1d(yl, gap=3), key=lambda t: abs(t[0] - best[0]))[0]
            c_r = min(C.clusters_1d(yr, gap=3), key=lambda t: abs(t[0] - best[0]))[0]
            slope = (c_r - c_l) / 12.0

        h_px = half_px + 0.5 + abs(slope) * max(xa.residual_px, 0.5)
        out.append(
            C.reading_point(
                x_val=float(xv),
                col=pc,
                row=row_abs,
                y_axis=ya,
                half_px=h_px,
                note=f"{label} row={row_idx}: stroke={2*half_px:.1f}px, res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
                clamp_min=0.0,
                clamp_max=float(max(y_tick_vals)),
            )
        )
    return out
