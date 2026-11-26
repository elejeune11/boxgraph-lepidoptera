"""Calibrated reader for Fig. 4B and Fig. 4D subdiffraction TIRFM normal intensity line scans."""

from __future__ import annotations

import numpy as np

import _calib as C


def _calibrate(image: np.ndarray):
    dark, gray, h_lines, v_lines = C.find_l_spines(
        image, thresh=100, row_frac=0.35, col_frac=0.35, row_range=(100, 1950), col_range=(250, 2150)
    )
    x_spine = max(h_lines, key=lambda t: t[0])
    y_spine = min(v_lines, key=lambda t: t[0])

    # Inward y ticks: 20, 15, 10, 5, 0 (x 10^2 a.u. -> 2000, 1500, 1000, 500, 0 a.u.)
    yt_raw = C.find_ticks(dark, y_spine, (180, int(x_spine[1]) + 5), "y", side=+1, band=(4, 24), frac=0.35, max_width=45)
    yt = C.match_ticks_regular(yt_raw, 5, tol=8.0, min_step=300.0)
    y_axis = C.verify_axis("linear", list(zip(yt, [2000.0, 1500.0, 1000.0, 500.0, 0.0])), max_residual_px=3.0)

    # Inward x ticks: 200, 400, 600, 800, 1000 nm
    xt_raw = C.find_ticks(dark, x_spine, (int(y_spine[2]) + 100, image.shape[1] - 50), "x", side=-1, band=(4, 24), frac=0.35, max_width=45)
    xt = C.match_ticks_regular(xt_raw, 5, tol=8.0, min_step=250.0)
    x_axis = C.verify_axis("linear", list(zip(xt, [200.0, 400.0, 600.0, 800.0, 1000.0])), max_residual_px=3.0)
    return x_axis, y_axis, x_spine, y_spine


def read_series(image: np.ndarray, channel: str, at: list[float]) -> list[dict]:
    x_axis, y_axis, x_spine, y_spine = _calibrate(image)
    R = image[:, :, 0].astype(np.float32)
    G = image[:, :, 1].astype(np.float32)
    B = image[:, :, 2].astype(np.float32)

    if channel == "og_actin":
        # Green Gaussian fit curve
        mask = (G > 75) & (G - R > 30) & (G - B > 15)
    elif channel == "phalloidin":
        # Magenta/pink Gaussian fit curve
        mask = (R > 150) & (B > 90) & (R - G > 45)
    else:
        raise ValueError(f"Unknown channel {channel!r}")

    row_min = 250
    row_max = int(x_spine[1]) - 5
    col_min = int(y_spine[2]) + 5
    mask[:row_min, :] = False
    mask[row_max + 1 :, :] = False
    mask[:, :col_min] = False

    return [
        C.sample_solid_curve_y(
            mask,
            x_axis,
            y_axis,
            float(xv),
            row_min,
            row_max,
            search_dx=8,
            note_prefix=f"Fig. 4 {channel}; ",
        )
        for xv in at
    ]
