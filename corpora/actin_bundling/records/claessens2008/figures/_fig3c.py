"""Shared calibration and curve/marker readers for Fig. 3C (p. 3): axial SAXS I(q_z) and four-sphere fits."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy import ndimage
from boxgraph.digitize import Axis

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import find_frame, find_ticks, sample_curve_at, value_interval


def calibrate_fig3c(image: np.ndarray):
    h, v = find_frame(image)
    r_top, r_bot = h[0][2], h[-1][2]
    c_left, c_right = v[0][2], v[-1][2]

    xt = find_ticks(image, r_bot, int(c_left) + 8, int(c_right) - 8, "x", -1, offset_lo=8, offset_hi=25)
    if len(xt) != 4:
        raise RuntimeError(f"Fig. 3C expected 4 x-ticks (1.1..1.4), found {len(xt)}")
    x_axis = Axis("linear", list(zip([p for p, _, _ in xt], [1.1, 1.2, 1.3, 1.4])))

    yt_l = find_ticks(image, c_left, int(r_top) + 8, int(r_bot) - 8, "y", +1, offset_lo=8, offset_hi=28)
    yt_r = find_ticks(image, c_right, int(r_top) + 8, int(r_bot) - 8, "y", -1, offset_lo=8, offset_hi=28)
    yt = [p for p, _, d in yt_r if d >= 26 and any(abs(p - pl) < 3.0 for pl, _, _ in yt_l)]
    if len(yt) != 7:
        raise RuntimeError(f"Fig. 3C expected 7 matched y-ticks (0.7..0.1), found {len(yt)}: {yt}")
    y_axis = Axis("linear", list(zip(yt, [0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1])))
    return x_axis, y_axis, (c_left, r_top, c_right, r_bot)


def _color_mask(image: np.ndarray, box: np.ndarray, series: str) -> np.ndarray:
    rgb = image.astype(float)
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    if series == "R0p05":
        return box & (R > 160.0) & (G < 80.0) & (B < 80.0)
    if series == "R0p1":
        return box & (B > 160.0) & (R < 80.0) & (G < 100.0)
    if series == "R0p5":
        return box & (G > 90.0) & (R < 70.0) & (B < 70.0)
    raise ValueError(f"unknown Fig. 3C series {series!r}")


def _fit_mask(mask: np.ndarray, radius: int = 4, min_comp: int = 120) -> np.ndarray:
    yy, xx = np.ogrid[-radius : radius + 1, -radius : radius + 1]
    dsk = (xx * xx + yy * yy) <= radius * radius
    opened = ndimage.binary_opening(mask, structure=dsk)
    lab, nc = ndimage.label(opened)
    clean = np.zeros_like(opened)
    for i in range(1, nc + 1):
        if (lab == i).sum() >= min_comp:
            clean |= lab == i
    return clean


def read_fig3c_fit(image: np.ndarray, at: list[float], series: str):
    x_axis, y_axis, (c_left, r_top, c_right, r_bot) = calibrate_fig3c(image)
    r0, r1 = int(r_top) + 8, int(r_bot) - 8
    c0, c1 = int(c_left) + 8, int(c_right) - 8
    box = np.zeros(image.shape[:2], bool)
    box[r0:r1, c0:c1] = True

    raw = _color_mask(image, box, series)
    fit = _fit_mask(raw, radius=4, min_comp=120)
    return [
        sample_curve_at(fit, x_axis, y_axis, qz, r0, r1, window=2, extra_px=1.5, note_prefix=f"fit {series}; ")
        for qz in at
    ]


def read_fig3c_data(image: np.ndarray, at: list[float], series: str):
    """Read measured SAXS markers at requested q_z positions by isolating marker envelopes
    (subtracting the thick continuous fit curve where separated, or centering the marker stroke envelope
    where markers and fit curve overlap)."""
    x_axis, y_axis, (c_left, r_top, c_right, r_bot) = calibrate_fig3c(image)
    r0, r1 = int(r_top) + 8, int(r_bot) - 8
    c0, c1 = int(c_left) + 8, int(c_right) - 8
    box = np.zeros(image.shape[:2], bool)
    box[r0:r1, c0:c1] = True

    raw = _color_mask(image, box, series)
    fit = _fit_mask(raw, radius=4, min_comp=120)
    # Dilate fit slightly so marker strokes sitting well away from the fit curve stand out cleanly
    fit_dil = ndimage.binary_dilation(fit, iterations=4)
    markers_only = raw & ~fit_dil

    cal_y = max(y_axis.residual_px, 0.5)
    readings = []
    for qz in at:
        c_exact = x_axis.to_pixel(qz)
        c_int = int(round(c_exact))
        # Examine a +/- 6 px column window (~half a marker width) to capture full marker height
        w_lo, w_hi = max(c0, c_int - 6), min(c1, c_int + 7)
        # First check if there is a clean marker-only cluster (separated from fit curve, e.g. wings/trough)
        m_patch = markers_only[r0:r1, w_lo:w_hi]
        f_patch = fit[r0:r1, w_lo:w_hi]
        r_patch = raw[r0:r1, w_lo:w_hi]

        m_rows = np.nonzero(m_patch.any(axis=1))[0]
        f_rows = np.nonzero(f_patch.any(axis=1))[0]
        r_rows = np.nonzero(r_patch.any(axis=1))[0]

        if len(r_rows) == 0:
            readings.append({"x": float(qz), "y": None, "note": f"data {series}: no pixels at q_z={qz:g}"})
            continue

        # Split marker-only rows into vertical clusters (gap > 22 px separates marker from fit edge artifact)
        chosen_row = None
        half_px = 6.0 + cal_y
        mode_note = ""
        if len(m_rows) > 0:
            splits = np.where(np.diff(m_rows) > 22)[0] + 1
            clusters = [g + r0 for g in np.split(m_rows, splits)]
            # Keep clusters that have substantial marker pixels (>= 20 px in the 13-col window)
            valid = []
            for g in clusters:
                area = int(m_patch[g[0] - r0 : g[-1] - r0 + 1, :].sum())
                if area >= 20:
                    valid.append((g, area))
            if valid:
                # If fit curve is well below (> 15 px lower in row, i.e. smaller I), pick the marker cluster above it;
                # otherwise pick the cluster containing the most marker pixels (largest hollow marker body)
                if len(f_rows) > 0:
                    f_med = float(np.median(f_rows) + r0)
                    above = [(g, a) for g, a in valid if float(g.mean()) < f_med - 15.0]
                    pick_g = max(above, key=lambda t: t[1])[0] if above else max(valid, key=lambda t: t[1])[0]
                else:
                    pick_g = max(valid, key=lambda t: t[1])[0]
                chosen_row = float(0.5 * (pick_g[0] + pick_g[-1]) + 0.5)
                span = float(pick_g[-1] - pick_g[0] + 1)
                # Note: if markers overlap the fit curve, also check full raw envelope around chosen_row
                if len(f_rows) > 0 and abs(chosen_row - (float(np.median(f_rows)) + r0)) < 25.0:
                    # Markers and fit coincide: center on full raw stroke envelope near chosen_row
                    near_r = r_rows[np.abs((r_rows + r0) - chosen_row) <= 22.0] + r0
                    chosen_row = float(0.5 * (near_r.min() + near_r.max()) + 0.5)
                    span = float(near_r.max() - near_r.min() + 1)
                    mode_note = "marker+fit envelope"
                else:
                    mode_note = "separated marker"
                half_px = max(4.0, 0.35 * min(span, 30.0)) + cal_y
        if chosen_row is None:
            chosen_row = float(0.5 * (r_rows.min() + r_rows.max()) + 0.5 + r0)
            half_px = 6.0 + cal_y
            mode_note = "raw envelope"

        y, lo, hi = value_interval(y_axis, chosen_row, half_px)
        readings.append({
            "x": float(qz),
            "y": y,
            "lo": lo,
            "hi": hi,
            "px": [float(c_exact), chosen_row],
            "px_lo": chosen_row - half_px,
            "px_hi": chosen_row + half_px,
            "note": f"data {series} ({mode_note}); cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f}) px",
        })
    return readings
