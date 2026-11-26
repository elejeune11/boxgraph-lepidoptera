"""Calibrated reader for Fig. S2F (1-Cumulative Frequency single-exponential fit curves at 50, 125, 500 nM fascin)."""

from __future__ import annotations

import numpy as np

import _calib as C


def _calibrate(image: np.ndarray):
    dark, gray, h_lines, v_lines = C.find_l_spines(
        image, thresh=100, row_frac=0.30, col_frac=0.30, row_range=(40, 2400), col_range=(200, image.shape[1] - 20)
    )
    x_spine = max(h_lines, key=lambda t: t[0])
    y_spine = min(v_lines, key=lambda t: t[0])

    # Inward y ticks: 1.0, [0.8 occluded by red curve], 0.6, 0.4, 0.2, 0.0
    yt_raw = C.find_ticks(dark, y_spine, (165, int(x_spine[1]) - 40), "y", side=+1, band=(4, 24), frac=0.35, max_width=45)
    # Take the 4 regular ticks [0.6, 0.4, 0.2, 0.0] in the lower part + 1.0 at ~198.5
    yt_low = C.match_ticks_regular([t for t in yt_raw if t > 850], 4, tol=8.0, min_step=340.0)
    step = (yt_low[-1] - yt_low[0]) / 3.0
    t_top = min(yt_raw, key=lambda t: abs(t - (yt_low[0] - 2.0 * step)))
    marks_y = [(t_top, 1.0), (yt_low[0], 0.6), (yt_low[1], 0.4), (yt_low[2], 0.2), (yt_low[3], 0.0)]
    y_axis = C.verify_axis("linear", marks_y, max_residual_px=4.0)

    # Inward x ticks every 0.5 um: 0.0 .. 9.5 um (20 ticks)
    xt_raw = C.find_ticks(dark, x_spine, (int(y_spine[0]) - 5, image.shape[1] - 20), "x", side=-1, band=(4, 24), frac=0.35, max_width=45)
    xt = C.match_ticks_regular(xt_raw, 20, tol=8.0, min_step=95.0)
    x_axis = C.verify_axis("linear", list(zip(xt, [0.5 * k for k in range(20)])), max_residual_px=4.5)
    return x_axis, y_axis, x_spine, y_spine


def _find_legend_box(image: np.ndarray, y_spine) -> tuple[int, int, int, int]:
    dark = image.astype(np.float32).mean(axis=2) < 100
    c_lo = int(y_spine[2]) + 300
    sub = dark[150:850, c_lo:]
    h_rows = np.where(sub.sum(axis=1) > 1000)[0] + 150
    if len(h_rows) < 2:
        return (0, 0, 0, 0)
    r0, r1 = int(h_rows.min()), int(h_rows.max())
    v_cols = np.where(dark[r0 : r1 + 1, c_lo:].sum(axis=0) > 0.7 * (r1 - r0))[0] + c_lo
    if len(v_cols) < 2:
        return (0, 0, 0, 0)
    return (int(v_cols.min()) - 8, r0 - 8, int(v_cols.max()) + 8, r1 + 8)


def read_series(image: np.ndarray, conc: str, at: list[float]) -> list[dict]:
    x_axis, y_axis, x_spine, y_spine = _calibrate(image)
    R = image[:, :, 0].astype(np.float32)
    G = image[:, :, 1].astype(np.float32)
    B = image[:, :, 2].astype(np.float32)

    if conc == "50nM":
        # Red curve
        mask = (R > 140) & (R - G > 55) & (R - B > 55)
    elif conc == "125nM":
        # Light sky-blue / cyan curve
        mask = (B > 170) & (G > 145) & (R < 155) & (B - R > 45)
    elif conc == "500nM":
        # Dark navy blue curve
        mask = (B > 90) & (B - R > 35) & (G < 130) & (R < 95)
    else:
        raise ValueError(f"Unknown conc {conc!r}")

    row_min = 120
    row_max = int(x_spine[1]) - 5
    col_min = int(y_spine[2]) + 5
    mask[:row_min, :] = False
    mask[row_max + 1 :, :] = False
    mask[:, :col_min] = False

    lc0, lr0, lc1, lr1 = _find_legend_box(image, y_spine)
    if lc1 > lc0 and lr1 > lr0:
        mask[lr0:lr1, lc0:lc1] = False

    return [
        C.sample_solid_curve_y(
            mask,
            x_axis,
            y_axis,
            float(xv),
            row_min,
            row_max,
            search_dx=10,
            note_prefix=f"Fig. S2F {conc}; ",
        )
        for xv in at
    ]
