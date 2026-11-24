"""Shared calibrated axis detection and curve reading helpers for SI Fig. S6 (PDF p. 20 / SI-11)."""
from __future__ import annotations

import math
import numpy as np
from boxgraph import digitize as D


def _clusters(mask: np.ndarray, min_len: int = 1, gap: int = 2) -> list[tuple[int, int, float]]:
    idx = np.where(mask)[0]
    if len(idx) == 0:
        return []
    splits = np.where(np.diff(idx) > gap)[0] + 1
    out = []
    for g in np.split(idx, splits):
        if len(g) >= min_len:
            out.append((int(g[0]), int(g[-1]), float(g.mean())))
    return out


def _color_mask(img: np.ndarray, t_key: str) -> np.ndarray:
    r = img[:, :, 0].astype(float)
    g = img[:, :, 1].astype(float)
    b = img[:, :, 2].astype(float)
    if t_key == "t10":
        return (r > 85.0) & (r < 175.0) & (g < 55.0) & (b < 55.0)
    if t_key == "t100":
        return (r > 210.0) & (g > 135.0) & (g < 230.0) & (b < 85.0)
    if t_key == "t1000":
        return (r < 140.0) & (g > 195.0) & (b > 140.0) & ((g - r) > 80.0)
    if t_key == "t10000":
        return (r < 85.0) & (g < 105.0) & (b > 205.0)
    if t_key == "t30000":
        return (r < 55.0) & (g < 55.0) & (b > 80.0) & (b < 175.0)
    raise ValueError(t_key)


def calibrate_s6(img: np.ndarray, panel: str):
    """Calibrate Fig. S6 panel ('a_left', 'a_right', 'b_left', 'b_right') from spines & inward ticks."""
    gray = img.astype(float).mean(axis=2)
    sat = img.max(axis=2).astype(float) - img.min(axis=2).astype(float)
    black = (gray < 140.0) & (sat < 30.0)
    vsp = _clusters(black.mean(axis=0) > 0.55)
    hsp = _clusters(black.mean(axis=1) > 0.55)
    c0, c1 = vsp[0][2], vsp[-1][2]
    r0, r1 = hsp[0][2], hsp[-1][2]
    c0i, c1i, r0i, r1i = int(round(c0)), int(round(c1)), int(round(r0)), int(round(r1))

    # Bottom inward x ticks: 0.5, 1.0, 1.5 (plus left spine 0.0 and right spine 2.0)
    xt_raw = _clusters(black[r1i - 42:r1i - 15, c0i + 150:c1i - 150].mean(axis=0) > 0.70)
    xt = [m + c0i + 150 for _, _, m in xt_raw]
    if len(xt) != 3:
        raise RuntimeError(f"Fig. S6 ({panel}): expected 3 interior x ticks (0.5, 1.0, 1.5), got {xt}")
    xa = D.Axis("linear", [(c0, 0.0), (xt[0], 0.5), (xt[1], 1.0), (xt[2], 1.5), (c1, 2.0)])

    # Left inward y ticks (upper ticks above the y=0 crowded baseline)
    yt_raw = _clusters(black[r0i + 100:r1i - 350, c0i + 15:c0i + 42].mean(axis=1) > 0.70)
    yt = [m + r0i + 100 for _, _, m in yt_raw]
    expected_y = {
        "a_left": [400.0, 300.0, 200.0, 100.0],
        "a_right": [15.0, 10.0, 5.0],
        "b_left": [600.0, 400.0, 200.0],
        "b_right": [10.0, 8.0, 6.0, 4.0, 2.0],
    }[panel]
    if len(yt) != len(expected_y):
        raise RuntimeError(f"Fig. S6 ({panel}): expected {len(expected_y)} upper y ticks, got {yt}")
    ya = D.Axis("linear", list(zip(yt, expected_y)))
    return xa, ya, (c0, r0, c1, r1)


def read_s6_curve(img: np.ndarray, at: list[float], panel: str, t_key: str) -> list[dict]:
    """Read one time-slice curve (t10, t100, t1000, t10000, t30000) in a Fig. S6 panel."""
    xa, ya, (c0, r0, c1, r1) = calibrate_s6(img, panel)
    bw, bh = c1 - c0, r1 - r0
    mask = _color_mask(img, t_key)

    clean = mask.copy()
    clean[:int(r0) + 10, :] = False
    clean[int(r1) - 10:, :] = False
    clean[:, :int(c0) + 6] = False
    clean[:, int(c1) - 6:] = False
    # Exclude legend in top-right corner
    clean[:int(r0 + 0.42 * bh), int(c0 + 0.54 * bw):] = False

    r_zero = ya.to_pixel(0.0)
    r_min_box = int(r0) + 10
    r_max_box = min(int(r1) - 10, int(round(r_zero + 65)))

    out = []
    for xv in at:
        c = xa.to_pixel(xv)
        # Try tight window first (+/- 10 px), widen to +/- 24 px only if occluded
        found = False
        for win in (10, 22, 36):
            ca = max(int(c0) + 6, int(round(c - win)))
            cb = min(int(c1) - 6, int(round(c + win)))
            sub = clean[r_min_box:r_max_box, ca:cb + 1]
            rr, cc = np.nonzero(sub)
            if len(rr) >= 12:
                rr_abs = rr + r_min_box
                # Look at vertical clusters in this narrow slice:
                # If a sharp spike has open circles & steep connecting lines, take median/midpoint of the primary band.
                # On steep flanks where circles span > 150 px vertically in a narrow slice, use the exact center
                # columns closest to c to avoid dragging the peak down!
                near_c = np.abs((cc + ca) - c) <= min(win, 8)
                if near_c.sum() >= 8:
                    r_use = rr_abs[near_c]
                else:
                    r_use = rr_abs
                # On needle-sharp spikes (vertical stem > 220 px in a narrow slice), a local peak column
                # has its open circle (outer radius ~45 px) at the very top (r_min), while the vertical stem
                # below would pull the raw pixel median down. Detect local peak columns purely from the image!
                v_span = float(r_use.max() - r_use.min())
                r_top_local = float(r_use.min())
                c_nbr_a = max(int(c0) + 6, int(round(c - 48)))
                c_nbr_b = min(int(c1) - 6, int(round(c + 48)))
                nbr_rr, nbr_cc = np.nonzero(clean[r_min_box:r_max_box, c_nbr_a:c_nbr_b + 1])
                is_local_peak = False
                if len(nbr_rr) > 0:
                    nbr_r_min = float(nbr_rr.min() + r_min_box)
                    cap_cols = nbr_cc[nbr_rr <= (nbr_rr.min() + 18)] + c_nbr_a
                    c_cap = float(np.median(cap_cols)) if len(cap_cols) else -1e9
                    if (r_top_local <= nbr_r_min + 12.0) and (abs(c - c_cap) <= 14.0):
                        is_local_peak = True
                if v_span > 220.0 and is_local_peak and panel.endswith("_left"):
                    r_center = r_top_local + 45.0
                    q25, q75 = r_center - 20.0, r_center + 20.0
                else:
                    q25, q50, q75 = np.percentile(r_use, [20, 50, 80])
                    r_center = float(q50)
                y_raw = ya.to_value(r_center)
                if y_raw < 0.0:
                    y_raw = 0.0
                    r_center = r_zero
                iqr_half = 0.5 * (q75 - q25)
                half_px = max(14.0, min(iqr_half, 48.0) + ya.residual_px + 4.0)
                if iqr_half > 80.0:
                    half_px = min(iqr_half * 0.75, 110.0)
                hi = ya.to_value(r_center - half_px)
                lo = max(0.0, ya.to_value(r_center + half_px))
                if y_raw < lo:
                    lo = y_raw
                out.append({
                    "x": float(xv), "y": float(y_raw), "lo": float(lo), "hi": float(hi),
                    "px": [float(c), float(r_center)],
                    "px_lo": float(r_center - half_px),
                    "px_hi": float(min(r_zero + 25.0, r_center + half_px)),
                    "note": f"x_res={xa.residual_px:.2f}px y_res={ya.residual_px:.2f}px win={win}px",
                })
                found = True
                break
        if not found:
            out.append({"x": float(xv), "y": None, "note": f"{t_key} occluded or absent at x={xv}"})
    return out
