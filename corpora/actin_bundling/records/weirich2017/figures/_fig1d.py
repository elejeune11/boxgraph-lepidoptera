"""Shared calibrated reader for Fig. 1D (actin/filamin phase diagram: filamin % vs filament length)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage, signal
import _calib as C


def calibrate(image: np.ndarray):
    dark, gray, hl, vl = C.find_box_spines(image, row_frac=0.5, col_frac=0.5)
    ht, hb = hl[0], hl[-1]
    vl0, vr0 = vl[0], vl[-1]
    xt = C.find_ticks(dark, hb, (int(vl0[2]) + 5, int(vr0[1]) - 5), "x", side=-1, band=(8, 45))
    yt = C.find_ticks(dark, vl0, (int(ht[2]) + 5, int(hb[1]) - 5), "y", side=+1, band=(8, 45))
    x_vals = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    if len(xt) != len(x_vals):
        raise RuntimeError(f"Fig1D: expected {len(x_vals)} x ticks, got {len(xt)}: {xt}")
    yt = C.match_ticks_regular(yt, 6, tol=4.0)
    xa = C.verify_axis("log10", [(vl0[0], 0.1)] + list(zip(xt, x_vals)), max_residual_px=1.5)
    ya = C.verify_axis("linear", list(zip(yt, [10.0, 8.0, 6.0, 4.0, 2.0, 0.0])), max_residual_px=1.5)
    return xa, ya, int(ht[2]) + 10, int(hb[1]) - 5, int(vl0[2]) + 10, int(vr0[1]) - 5


def _matched_peaks(mask: np.ndarray, kernel: np.ndarray, thresh_frac: float, min_dist: int = 45) -> list[tuple[float, float, float]]:
    """Find local maxima of normalized cross-correlation with a binary marker kernel."""
    # Downsample 2x for speed (still 600 dpi -> 0.5 px precision after upscaling!)
    m2 = mask[::2, ::2].astype(np.float32)
    k2 = kernel[::2, ::2].astype(np.float32)
    k_sum = float(k2.sum())
    # Penalize background slightly so solid blobs don't falsely match hollow shapes
    k_norm = k2 / k_sum
    score = signal.fftconvolve(m2, k_norm[::-1, ::-1], mode="same")
    # Local maximum filter
    sz = max(3, min_dist // 2)
    loc_max = score == ndimage.maximum_filter(score, size=(sz, sz))
    cand = loc_max & (score >= thresh_frac)
    rr, cc = np.nonzero(cand)
    peaks = sorted([(float(score[r, c]), 2.0 * c + 0.5, 2.0 * r + 0.5) for r, c in zip(rr, cc)], reverse=True)
    # Greedy non-maximum suppression at min_dist in full-res px
    kept: list[tuple[float, float, float]] = []
    for sc, c_f, r_f in peaks:
        if all((c_f - kc) ** 2 + (r_f - kr) ** 2 >= min_dist**2 for _, kc, kr in kept):
            kept.append((sc, c_f, r_f))
    return kept


def _format_points(pts: list[tuple[float, float]], xa, ya, label: str, half_px: float = 4.0) -> list[dict]:
    out = []
    hp = half_px + ya.residual_px
    for c_f, r_f in pts:
        xv = float(xa.to_value(c_f))
        yv = float(ya.to_value(r_f))
        ya_, yb_ = float(ya.to_value(r_f - hp)), float(ya.to_value(r_f + hp))
        out.append(
            {
                "x": xv,
                "y": yv,
                "lo": min(ya_, yb_),
                "hi": max(ya_, yb_),
                "px": [float(c_f), float(r_f)],
                "px_lo": float(r_f - hp),
                "px_hi": float(r_f + hp),
                "note": f"Fig1D {label}, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
            }
        )
    out.sort(key=lambda d: (round(d["x"], 4), d["y"]))
    return out


def read_class(image: np.ndarray, which: str) -> list[dict]:
    xa, ya, r0, r1, c0, c1 = calibrate(image)
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    inb = np.zeros(image.shape[:2], bool)
    inb[r0:r1, c0:c1] = True
    # Exclude text labels 'LIQUID droplets' and 'VISCOELASTIC networks'
    inb[int(ya.to_pixel(9.0)) : int(ya.to_pixel(6.5)), int(xa.to_pixel(0.14)) : int(xa.to_pixel(0.55))] = False
    inb[int(ya.to_pixel(9.8)) : int(ya.to_pixel(7.6)), int(xa.to_pixel(1.35)) : int(xa.to_pixel(9.2))] = False

    if which == "liquid_droplets":
        # Light blue / cyan filled circles (diameter ~116 px)
        cyan = inb & (b - r > 80) & (g - r > 50) & (b > 170)
        # Distance transform >= 43 separates closely spaced touching circles cleanly into individual cores
        cores = ndimage.distance_transform_edt(cyan) >= 43.0
        pts = [(bl["c"], bl["r"]) for bl in C.find_blobs(cores, min_area=120)]
        return _format_points(pts, xa, ya, "liquid droplet", half_px=4.0)

    if which == "fluid_bundles":
        # Dark blue open rings (RGB ~(41, 41, 118), outer diam ~188 px, stroke ~24 px)
        dblue = inb & (b - r > 35) & (b - g > 35) & (g < 95) & (r < 95) & (b > 80)
        yy, xx = np.ogrid[-100:101, -100:101]
        rad = np.sqrt(xx * xx + yy * yy)
        ring = (rad >= 71.0) & (rad <= 95.0)
        peaks = _matched_peaks(dblue, ring, thresh_frac=0.45, min_dist=80)
        pts = [(c_f, r_f) for _, c_f, r_f in peaks]
        return _format_points(pts, xa, ya, "fluid bundle", half_px=5.0)

    if which == "viscoelastic_networks":
        # Medium gray X crosses (RGB ~(115, 115, 115), ~130x130 px, arm width ~24 px)
        inb_x = inb.copy()
        inb_x[r0:r1, c0 : image.shape[1] - 5] = True
        inb_x[int(ya.to_pixel(9.8)) : int(ya.to_pixel(7.6)), int(xa.to_pixel(1.35)) : int(xa.to_pixel(9.2))] = False
        medgray = inb_x & (r > 75) & (r < 170) & (np.abs(r - g) < 16) & (np.abs(r - b) < 16)
        yy, xx = np.ogrid[-70:71, -70:71]
        x_kern = ((np.abs(xx - yy) <= 14) | (np.abs(xx + yy) <= 14)) & (np.abs(xx) <= 60) & (np.abs(yy) <= 60)
        peaks = _matched_peaks(medgray, x_kern, thresh_frac=0.60, min_dist=85)
        pts = [(c_f, r_f) for _, c_f, r_f in peaks]
        return _format_points(pts, xa, ya, "viscoelastic network", half_px=5.0)

    if which == "gas":
        # Black '+' symbols in bottom-left gas region (x < 0.68, y < 1.8)
        gas_box = np.zeros(image.shape[:2], bool)
        gas_box[int(ya.to_pixel(1.8)) : r1, c0 : int(xa.to_pixel(0.68))] = True
        blk = gas_box & (r < 45) & (g < 45) & (b < 45)
        yy, xx = np.ogrid[-82:83, -82:83]
        plus_kern = ((np.abs(xx) <= 14) & (np.abs(yy) <= 74)) | ((np.abs(yy) <= 14) & (np.abs(xx) <= 74))
        peaks = _matched_peaks(blk, plus_kern, thresh_frac=0.68, min_dist=70)
        pts = [(c_f, r_f) for _, c_f, r_f in peaks]
        return _format_points(pts, xa, ya, "gas single filaments", half_px=5.0)

    raise ValueError(which)
