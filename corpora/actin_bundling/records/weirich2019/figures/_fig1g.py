"""Shared calibrated reader for Main Fig. 1G Upper (high % active) and Lower (lower % active)."""

from __future__ import annotations

import numpy as np
from _calib import (
    find_blobs,
    find_box_spines,
    find_ticks,
    make_point_reading,
    match_ticks_regular,
    verify_axis,
)


def _calibrate(image: np.ndarray):
    dark, gray, hs, vs = find_box_spines(image, thresh=130, row_frac=0.30, col_frac=0.45)
    hs = [s for s in hs if (s[2] - s[1]) < 25]
    vs = [s for s in vs if (s[2] - s[1]) < 25]
    top, bot = hs[0], hs[-1]
    left, right = vs[0], vs[-1]

    # Top spine always has 6 clean unobstructed x-ticks at 2, 4, 6, 8, 10, 12 nM; left spine is 0 nM.
    xt_raw = find_ticks(dark, top, (left[2] + 10, right[1] - 10), "x", side=+1, band=(5, 25), frac=0.5)
    xt = match_ticks_regular(xt_raw, 6, tol=5.0)
    x_marks = [(left[0], 0.0)] + [(xt[i], float(2 * (i + 1))) for i in range(6)]
    x_axis = verify_axis("linear", x_marks, max_residual_px=2.0)

    # Left and right spines have 4 ticks at 0.8, 0.6, 0.4, 0.2; top = 1.0, bot = 0.0.
    yt_raw = find_ticks(dark, left, (top[2] + 10, bot[1] - 10), "y", side=+1, band=(5, 25), frac=0.5)
    yt = match_ticks_regular(yt_raw, 4, tol=5.0)
    y_marks = [(top[0], 1.0)] + [(yt[i], 0.8 - 0.2 * i) for i in range(4)] + [(bot[0], 0.0)]
    y_axis = verify_axis("linear", y_marks, max_residual_px=2.0)

    return x_axis, y_axis, (left, top, right, bot)


def _color_masks(image: np.ndarray, panel: str, spines):
    left, top, right, bot = spines
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)

    red = (R > 160) & (G < 80) & (B < 80)
    blue = (B > 100) & ((B - R) > 30) & ((B - G) > 30)
    purple = (R > 80) & (B > 80) & (G < 100) & (np.abs(R - B) < 50) & ((R - G) > 30)

    H, W = image.shape[:2]
    roi = np.zeros((H, W), dtype=bool)
    roi[max(0, top[1] - 40) : min(H, bot[2] + 65), left[2] + 10 : right[1] - 10] = True

    # Mask out legend box only (in Lower, legend sits above bottom spine so y=0 markers stay visible)
    if panel == "upper":
        roi[: top[2] + 400, right[1] - 700 :] = False
    else:
        roi[bot[1] - 420 : bot[1] - 60, right[1] - 700 :] = False

    return red & roi, blue & roi, purple & roi


def _merge_marker_clusters(mask: np.ndarray, col_tol: float = 55.0, min_area: int = 350):
    """Find markers and merge pieces split by an overlapping marker of another color."""
    blobs = find_blobs(mask, min_area=min_area)
    # Filter out small dotted connector dashes (markers have bounding span >= 35 px or area >= 700)
    blobs = [b for b in blobs if b["area"] >= 700 or max(b["w"], b["h"]) >= 40]
    blobs.sort(key=lambda b: b["c"])

    groups: list[list[dict]] = []
    for b in blobs:
        if groups and abs(b["c"] - np.mean([g["c"] for g in groups[-1]])) <= col_tol:
            groups[-1].append(b)
        else:
            groups.append([b])

    centers = []
    for g in groups:
        c0 = min(b["c0"] for b in g)
        c1 = max(b["c1"] for b in g)
        r0 = min(b["r0"] for b in g)
        r1 = max(b["r1"] for b in g)
        # Check if a single clean unobstructed marker vs partially overlapped
        if len(g) == 1 and abs(g[0]["w"] - g[0]["h"]) <= 12:
            # Clean symmetric marker: centroid is sub-pixel accurate
            centers.append((g[0]["c"], g[0]["r"], 4.0, f"blob {g[0]['w']}x{g[0]['h']}px"))
        else:
            # Partially occluded by overlapping marker: outer bounding-box center
            c_mid = 0.5 * (c0 + c1) + 0.5
            r_mid = 0.5 * (r0 + r1) + 0.5
            centers.append((c_mid, r_mid, 6.0, f"merged bbox {c1 - c0 + 1}x{r1 - r0 + 1}px (n={len(g)})"))
    return centers


def read_series(image: np.ndarray, panel: str, series: str) -> list[dict]:
    x_axis, y_axis, spines = _calibrate(image)
    red, blue, purple = _color_masks(image, panel, spines)

    if series == "centered":
        raw_centers = _merge_marker_clusters(red, col_tol=60.0, min_area=500)
    elif series == "bisected":
        raw_centers = _merge_marker_clusters(blue, col_tol=60.0, min_area=500)
    elif series == "separated":
        raw_centers = _merge_marker_clusters(purple, col_tol=60.0, min_area=500)
    else:
        raise ValueError(f"Unknown series {series!r}")

    # For Fig. 1G lower, if a blue 'x' at y=0 has its top occluded or only the central portion visible
    # inside the purple square, check its symmetry relative to the square center; all markers in Fig. 1G
    # have 3 concentrations per panel.
    if len(raw_centers) != 3:
        raise RuntimeError(f"Fig. 1G {panel}/{series}: expected 3 markers, found {len(raw_centers)}: {raw_centers}")

    cal_note = f"cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px"
    out = []
    for c, r, h_px, info in raw_centers:
        out.append(make_point_reading(c, r, x_axis, y_axis, half_px_y=h_px, note=f"{info}, {cal_note}"))
    return out
