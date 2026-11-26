"""Shared calibration and curve extraction for Fig. 1c (p. 2 / printed p. 790)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import _common as C


def calibrate_fig1c(image: np.ndarray):
    box = C.find_plot_box(image, dark_thr=120.0, sat_thr=35.0, min_h_frac=0.55, min_v_frac=0.55)
    dark = C.dark_neutral_mask(image, dark_thr=120.0, sat_thr=35.0)
    top, bot, left, right = box["top"], box["bot"], box["left"], box["right"]

    # Inward x ticks on bottom spine: 0.2, 0.4, 0.6, 0.8, 1.0
    # Note: Nature's art redraw nudged the outermost tick '1.0' ~32.5 px inward and y ticks '2.5'/'0.5'
    # ~24/8 px inward relative to the strictly uniform interior ticks (0.2..0.8 spaced 339.5 +/- 0.5 px;
    # 1.0..2.0 spaced 354.25 +/- 0.25 px). We calibrate from the uniform interior ticks and propagate
    # the outer-tick typesetting discrepancy honestly into the reading interval.
    xt_all = C.scan_ticks(dark, "x", bot[0] - 26, bot[0] - 6, left[1] + 5, right[0] - 5, min_frac=0.6)
    yt_all = C.scan_ticks(dark, "y", top[1] + 5, bot[0] - 5, left[1] + 6, left[1] + 26, min_frac=0.6)

    xt_reg = C.select_regular_ticks(xt_all, 4, tol=3.0)
    xa = C.make_axis("linear", list(zip(xt_reg, [0.2, 0.4, 0.6, 0.8])), max_res_px=2.0)

    yt_reg = C.select_regular_ticks(yt_all, 3, tol=3.0)
    ya = C.make_axis("linear", list(zip(yt_reg, [2.0, 1.5, 1.0])), max_res_px=2.0)

    # Measure max deviation of outermost printed ticks from uniform interior scale for honest uncertainty
    extra_y_px = 0.0
    if len(yt_all) >= 5:
        y_top_dev = abs(min(yt_all) - ya.to_pixel(2.5))
        y_bot_dev = abs(max(yt_all) - ya.to_pixel(0.5))
        extra_y_px = 0.5 * max(y_top_dev, y_bot_dev)

    return box, xa, ya, extra_y_px


def read_series(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    box, xa, ya, extra_y_px = calibrate_fig1c(image)
    top, bot, left, right = box["top"], box["bot"], box["left"], box["right"]
    rgb = image.astype(np.float32)
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]

    interior = np.zeros(image.shape[:2], dtype=bool)
    interior[top[1] + 6 : bot[0] - 6, left[1] + 6 : right[0] - 6] = True
    # Mask out upper-left legend block (Actin / Arp2/3 / CP) using calibrated coordinates
    leg_r1 = int(round(ya.to_pixel(1.38)))
    leg_c1 = int(round(xa.to_pixel(0.50)))
    interior[:leg_r1, :leg_c1] = False

    if series == "actin":
        mask = interior & (R > 165.0) & ((R - G) > 80.0) & ((R - B) > 80.0)
        lbl = "Fig. 1c Actin (red)"
    elif series == "cp":
        mask = interior & (B > 110.0) & ((B - R) > 45.0) & ((B - G) > 35.0)
        lbl = "Fig. 1c CP/Actin (blue)"
    elif series == "arp23":
        # Neutral black curve (exclude red and blue strokes and ticks)
        mask = interior & (rgb.max(axis=2) < 85.0)
        # Exclude inward bottom ticks and left ticks
        mask[bot[0] - 32 :, :] = False
        mask[:, : left[1] + 32] = False
        lbl = "Fig. 1c Arp2/3/Actin (black)"
    else:
        raise ValueError(f"Unknown series {series!r}")

    return C.trace_curve_at(
        mask=mask,
        x_axis=xa,
        y_axis=ya,
        at=at,
        r_min=top[1] + 6,
        r_max=bot[0] - 6,
        half_col=3,
        search_window_px=10,
        min_run=6,
        extra_y_px=extra_y_px,
        clamp_min=0.0,
        label=lbl,
    )
