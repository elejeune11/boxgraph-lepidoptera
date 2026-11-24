"""Shared calibrated reader for Fig. 3B (tactoid coalescence length vs time)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage, signal
import _calib as C


def calibrate(image: np.ndarray):
    dark, gray, hl, vl = C.find_box_spines(image, row_frac=0.5, col_frac=0.5)
    ht, hb = hl[0], hl[-1]
    vl0, vr0 = vl[0], vl[-1]
    xt = C.find_ticks(dark, hb, (int(vl0[2]) + 5, int(vr0[1]) - 5), "x", side=-1, band=(8, 45))
    xt = C.match_ticks_regular(xt, 6, tol=4.0)
    yt = C.find_ticks(dark, vr0, (int(ht[2]) + 5, int(hb[1]) - 5), "y", side=-1, band=(8, 45))
    yt = C.match_ticks_regular(yt, 6, tol=4.0)

    xa = C.verify_axis(
        "linear", [(vl0[0], -100.0)] + list(zip(xt, [-50.0, 0.0, 50.0, 100.0, 150.0, 200.0])) + [(vr0[0], 250.0)], 1.5
    )
    ya = C.verify_axis(
        "linear", [(ht[0], 10.0)] + list(zip(yt, [9.0, 8.0, 7.0, 6.0, 5.0, 4.0])) + [(hb[0], 3.0)], 1.5
    )
    return xa, ya, int(ht[2]) + 15, int(hb[1]) - 15, int(vl0[2]) + 5, int(vr0[1]) - 5


def _blue_dots(image: np.ndarray, xa, ya, r0: int, r1: int, c0: int, c1: int) -> list[tuple[float, float]]:
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    inb = np.zeros(image.shape[:2], bool)
    inb[r0:r1, c0:c1] = True
    blue = inb & (b - r > 80.0) & (b > 140.0)
    # Ignore black fit dashes / axis spines in background penalty so dashes crossing blue dots don't bias centers
    blk = (r < 60.0) & (g < 60.0) & (b < 60.0)
    bg = inb & (~blue) & (~blk)
    # Measured isolated blue dot radius is 28.0 px; greedy disk matching pursuit separates all 5-s overlapping frames
    yy, xx = np.ogrid[-28:29, -28:29]
    disk = (xx * xx + yy * yy) <= 27.5**2
    d2 = disk[::2, ::2].astype(np.float32)
    b2 = blue[::2, ::2].astype(np.float32)
    bg2 = bg[::2, ::2].astype(np.float32)
    pen = signal.fftconvolve(bg2, d2[::-1, ::-1], mode="same")
    uncov = b2.copy()
    d_sum = float(d2.sum())
    h2 = d2.shape[0] // 2
    d_mask = d2 > 0.5
    kept: list[tuple[float, float]] = []
    # Restrict centers inside [c0+14, c1-14] so spine-clipped edge fragments aren't shifted
    c_min_2, c_max_2 = (c0 + 14) // 2, (c1 - 14) // 2
    for _ in range(90):
        gain = signal.fftconvolve(uncov, d2[::-1, ::-1], mode="same")
        valid = (pen <= 0.14 * d_sum) & (gain >= 0.40 * d_sum)
        valid[:, :c_min_2] = False
        valid[:, c_max_2:] = False
        if not np.any(valid):
            break
        sc = np.where(valid, gain - 2.0 * pen, -1e9)
        ri, ci = np.unravel_index(int(np.argmax(sc)), sc.shape)
        if sc[ri, ci] < 0:
            break
        kept.append((2.0 * ci + 0.5, 2.0 * ri + 0.5))
        r_s0, c_s0 = ri - h2, ci - h2
        uncov[r_s0 : r_s0 + d2.shape[0], c_s0 : c_s0 + d2.shape[1]] *= ~d_mask
    return kept


def read_series(image: np.ndarray, at: list[float], which: str) -> list[dict]:
    xa, ya, r0, r1, c0, c1 = calibrate(image)

    if which in ("pre_upper", "pre_lower", "coalesced_data"):
        pts = _blue_dots(image, xa, ya, r0, r1, c0, c1)
        out = []
        half_px = 3.0 + ya.residual_px
        for c_m, r_m in pts:
            xv = float(xa.to_value(c_m))
            yv = float(ya.to_value(r_m))
            if which == "pre_upper" and not (xv <= 3.0 and yv > 4.65):
                continue
            if which == "pre_lower" and not (xv <= 3.0 and yv <= 4.65):
                continue
            if which == "coalesced_data" and not (xv > 3.0):
                continue
            ya_, yb_ = float(ya.to_value(r_m - half_px)), float(ya.to_value(r_m + half_px))
            out.append(
                {
                    "x": xv,
                    "y": yv,
                    "lo": min(ya_, yb_),
                    "hi": max(ya_, yb_),
                    "px": [c_m, r_m],
                    "px_lo": float(r_m - half_px),
                    "px_hi": float(r_m + half_px),
                    "note": f"Fig3B {which} dot, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
                }
            )
        out.sort(key=lambda d: d["x"])
        return out

    if which == "exp_fit":
        r = image[:, :, 0].astype(float)
        g = image[:, :, 1].astype(float)
        b = image[:, :, 2].astype(float)
        # Pure neutral black dashes for t > 10 s
        blk = (r < 50) & (g < 50) & (b < 50) & (np.abs(r - g) < 15) & (np.abs(r - b) < 15)
        m = np.zeros_like(blk)
        r_lo = int(ya.to_pixel(9.35))
        r_hi = int(ya.to_pixel(5.65))
        c_lo = int(xa.to_pixel(9.0))
        c_hi = int(xa.to_pixel(249.0))
        m[r_lo:r_hi, c_lo:c_hi] = blk[r_lo:r_hi, c_lo:c_hi]
        # Mask Li arrow + scale bar (t in 25..135, y > 8.6) and right-spine tick at y=6 (t > 242)
        m[: int(ya.to_pixel(8.55)), int(xa.to_pixel(25.0)) : int(xa.to_pixel(135.0))] = False
        m[:, int(xa.to_pixel(244.0)) :] = False
        return [
            C.sample_dashed_curve_y(
                m, xa, ya, float(x), r_lo, r_hi, max_search_dx=110, sigma_dx=45.0, deg=2, note_prefix="Fig3B exp fit: "
            )
            for x in at
        ]

    raise ValueError(which)
