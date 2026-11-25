"""Shared reader for Fig. 3(d) (PDF p. 4 / 128402-4): Cahn-Hilliard domain size a(t) vs. t (Shear & Rotational)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage

from _common import build_axis, find_box, match_ticks_both_sides, reading_from_px, scan_ticks


def calibrate_fig3d(image: np.ndarray):
    box = find_box(image)
    c0, r0, c1, r1 = box

    xt = [
        p
        for p, _ in match_ticks_both_sides(
            scan_ticks(image, box, "bottom_in", offset_min=5, offset_max=25, min_len=5, gray_thresh=200),
            scan_ticks(image, box, "top_in", offset_min=5, offset_max=25, min_len=5, gray_thresh=200),
        )
        if c0 + 50 < p < c1 - 50
    ]
    xa = build_axis("log10", [c0] + xt + [c1], [20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0], max_residual_px=2.0)

    yt = [
        p
        for p, _ in scan_ticks(image, box, "left_in", offset_min=5, offset_max=25, min_len=5, gray_thresh=200)
        if r0 + 50 < p < r1 - 50
    ]
    ya = build_axis("log10", yt + [r1], [20.0, 18.0, 16.0, 14.0, 12.0, 10.0], max_residual_px=2.0)
    return box, xa, ya


def _clean_curves_3d(image: np.ndarray, box: tuple[float, float, float, float]) -> np.ndarray:
    """Keep the connected curve+marker bands spanning the plot, erasing disconnected text ('Shear', 'Rotational', 'alpha = ...')."""
    c0, r0, c1, r1 = box
    dark = image.astype(np.float64).mean(axis=2) < 140.0
    r_lo, r_hi = int(r0) + 25, int(r1) - 25
    c_lo, c_hi = int(c0) + 12, int(c1) - 12
    interior = np.zeros_like(dark)
    interior[r_lo:r_hi, c_lo:c_hi] = dark[r_lo:r_hi, c_lo:c_hi]

    # Dilate slightly to connect adjacent markers to their fit line, keep only components spanning > 1000 px wide
    dil = ndimage.binary_dilation(interior, iterations=6)
    lbl, _ = ndimage.label(dil)
    slices = ndimage.find_objects(lbl)
    keep_mask = np.zeros_like(interior)
    for i, sl in enumerate(slices, start=1):
        if sl is None:
            continue
        if (sl[1].stop - sl[1].start) >= 1000:
            keep_mask[sl] |= lbl[sl] == i
    return interior & keep_mask


def read_series(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    box, xa, ya = calibrate_fig3d(image)
    c0, r0, c1, r1 = box
    clean = _clean_curves_3d(image, box)

    out = []
    half_w = 12
    r_lo, r_hi = int(r0) + 25, int(r1) - 25
    for t in at:
        col = xa.to_pixel(t)
        ci = int(round(col))
        strip = clean[r_lo:r_hi, max(int(c0) + 12, ci - half_w) : min(int(c1) - 12, ci + half_w + 1)]
        rr, _ = np.nonzero(strip)
        if len(rr) == 0:
            out.append({"x": float(t), "y": None, "found": False, "note": f"no pixels at t={t}"})
            continue
        # Split vertical clusters separated by >= 55 px white gap (upper = Shear, lower = Rotational)
        u_rows = np.sort(np.unique(rr))
        gaps = np.where(np.diff(u_rows) >= 55)[0]
        if len(gaps) == 0:
            out.append({"x": float(t), "y": None, "found": False, "note": f"curves not separated at t={t}"})
            continue
        # Largest vertical gap separates Shear (above) from Rotational (below)
        g_idx = gaps[int(np.argmax(np.diff(u_rows)[gaps]))]
        split_r = 0.5 * (u_rows[g_idx] + u_rows[g_idx + 1])
        if series == "shear":
            sel = rr < split_r
        elif series == "rotational":
            sel = rr > split_r
        else:
            raise ValueError(f"unknown series {series!r}")
        rows_abs = rr[sel] + r_lo + 0.5
        row_med = float(np.median(rows_abs))
        p10, p90 = float(np.percentile(rows_abs, 10)), float(np.percentile(rows_abs, 90))
        half_px = max(0.5 * (p90 - p10), 12.0)
        out.append(
            reading_from_px(
                xa,
                ya,
                col,
                row_med,
                half_px,
                x_val=float(t),
                note=f"Fig. 3(d) {series} median of {len(rows_abs)} px (half-span {half_px:.1f} px, cal res ({xa.residual_px:.2f}, {ya.residual_px:.2f}) px)",
            )
        )
    return out
