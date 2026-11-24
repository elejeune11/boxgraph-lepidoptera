"""Shared calibrated reader for Fig. 3C inset (relaxation time tau vs deformation length L_i - L_f)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage, signal
import _calib as C


def calibrate(image: np.ndarray):
    dark, gray, hl, vl = C.find_box_spines(image, row_frac=0.5, col_frac=0.5)
    ht, hb = hl[0], hl[-1]
    vl0, vr0 = vl[0], vl[-1]
    xt = C.find_ticks(dark, hb, (int(vl0[2]) + 5, int(vr0[1]) - 5), "x", side=-1, band=(6, 40))
    xt = C.match_ticks_regular(xt, 7, tol=4.0)
    yt = C.find_ticks(dark, vl0, (int(ht[2]) + 5, int(hb[1]) - 5), "y", side=+1, band=(6, 40))
    # Note: y=10 tick is shifted ~4 px upward in the published vector graphic relative to the regular 202-px grid
    if len(yt) != 6:
        raise RuntimeError(f"Fig3C inset: expected 6 interior y ticks (60..10), got {yt}")

    xa = C.verify_axis(
        "linear", [(vl0[0], 0.5)] + list(zip(xt, [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0])) + [(vr0[0], 4.5)], 1.5
    )
    ya = C.verify_axis("linear", list(zip(yt, [60.0, 50.0, 40.0, 30.0, 20.0, 10.0])) + [(hb[0], 0.0)], 3.5)
    return xa, ya, int(ht[2]) + 10, int(hb[1]) - 10, int(vl0[2]) + 10, int(vr0[1]) - 10


def _mask_legend(mask: np.ndarray, xa, ya) -> np.ndarray:
    out = mask.copy()
    # Legend box sits in top-left (x < 1.9, y > 38)
    out[: int(ya.to_pixel(38.0)), : int(xa.to_pixel(1.9))] = False
    return out


def read_series(image: np.ndarray, at: list[float], which: str) -> list[dict]:
    xa, ya, r0, r1, c0, c1 = calibrate(image)
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    inb = np.zeros(image.shape[:2], bool)
    inb[r0:r1, c0:c1] = True
    inb = _mask_legend(inb, xa, ya)

    if which == "10pct":
        blue = inb & (b - r > 50) & (b > 150)
        cores = ndimage.distance_transform_edt(blue) >= 20.0
        out = []
        half_px = 3.0 + ya.residual_px
        for bl in C.find_blobs(cores, min_area=100):
            c_m, r_m = bl["c"], bl["r"]
            yv = float(ya.to_value(r_m))
            ya_, yb_ = float(ya.to_value(r_m - half_px)), float(ya.to_value(r_m + half_px))
            out.append(
                {
                    "x": float(xa.to_value(c_m)),
                    "y": yv,
                    "lo": min(ya_, yb_),
                    "hi": max(ya_, yb_),
                    "px": [c_m, r_m],
                    "px_lo": float(r_m - half_px),
                    "px_hi": float(r_m + half_px),
                    "note": f"Fig3C inset 10% blue circle, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
                }
            )
        out.sort(key=lambda d: d["x"])
        return out

    if which == "5pct":
        red = inb & (r - g > 80) & (r - b > 80)
        # Matched filter with X kernel (~70x70 px, stroke ~10 px) separates slightly touching red crosses cleanly
        yy, xx = np.ogrid[-36:37, -36:37]
        xk = ((np.abs(xx - yy) <= 6) | (np.abs(xx + yy) <= 6)) & (np.abs(xx) <= 31) & (np.abs(yy) <= 31)
        m2 = red.astype(np.float32)
        k_norm = xk.astype(np.float32) / float(xk.sum())
        score = signal.fftconvolve(m2, k_norm[::-1, ::-1], mode="same")
        loc_max = (score == ndimage.maximum_filter(score, size=(35, 35))) & (score >= 0.45)
        rr, cc = np.nonzero(loc_max)
        peaks = sorted([(float(score[ri, ci]), float(ci) + 0.5, float(ri) + 0.5) for ri, ci in zip(rr, cc)], reverse=True)
        kept = []
        for sc, cf, rf in peaks:
            if all((cf - kc) ** 2 + (rf - kr) ** 2 >= 36.0**2 for _, kc, kr in kept):
                kept.append((sc, cf, rf))
        out = []
        half_px = 3.0 + ya.residual_px
        for _, c_m, r_m in kept:
            yv = float(ya.to_value(r_m))
            ya_, yb_ = float(ya.to_value(r_m - half_px)), float(ya.to_value(r_m + half_px))
            out.append(
                {
                    "x": float(xa.to_value(c_m)),
                    "y": yv,
                    "lo": min(ya_, yb_),
                    "hi": max(ya_, yb_),
                    "px": [c_m, r_m],
                    "px_lo": float(r_m - half_px),
                    "px_hi": float(r_m + half_px),
                    "note": f"Fig3C inset 5% red cross, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
                }
            )
        out.sort(key=lambda d: d["x"])
        return out

    if which == "linear_fit":
        # Pure black dashes, inset past inward ticks (45 px) and outside legend box
        blk = (r < 45) & (g < 45) & (b < 45) & (np.abs(r - g) < 15) & (np.abs(r - b) < 15)
        m = np.zeros_like(blk)
        m[r0 + 38 : r1 - 38, c0 + 38 : c1 - 38] = blk[r0 + 38 : r1 - 38, c0 + 38 : c1 - 38]
        m = _mask_legend(m, xa, ya)
        return [
            C.sample_dashed_curve_y(
                m, xa, ya, float(x), r0 + 38, r1 - 38, max_search_dx=180, sigma_dx=80.0, deg=1, note_prefix="Fig3C inset linear fit: "
            )
            for x in at
        ]

    raise ValueError(which)
