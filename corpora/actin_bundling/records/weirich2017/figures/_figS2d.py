"""Shared calibrated reader for SI Fig. S2D (aspect ratio L/r vs omega, log-log)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage
import _calib as C


def calibrate(image: np.ndarray):
    dark, gray, hl, vl = C.find_box_spines(image, row_frac=0.5, col_frac=0.5)
    ht, hb = hl[0], hl[-1]
    vl0, vr0 = vl[0], vl[-1]
    xt_major = C.find_ticks(dark, hb, (int(vl0[2]) + 10, int(vr0[1]) - 10), "x", side=-1, band=(8, 45))
    if len(xt_major) != 3:
        raise RuntimeError(f"FigS2D: expected 3 interior major x ticks (0.1, 0.5, 1.0), got {xt_major}")
    # Use right spine for major y ticks (4, 3, 2, 1) to avoid blue sphere icon touching left spine
    yt_major = C.find_ticks(dark, vr0, (int(ht[2]) + 10, int(hb[1]) - 10), "y", side=-1, band=(8, 45))
    if len(yt_major) != 4:
        raise RuntimeError(f"FigS2D: expected 4 major y ticks on right spine (4, 3, 2, 1), got {yt_major}")

    x_marks = list(zip(xt_major, [0.1, 0.5, 1.0])) + [(vr0[0], 5.0)]
    y_marks = list(zip(yt_major, [4.0, 3.0, 2.0, 1.0]))
    xa = C.verify_axis("log10", x_marks, max_residual_px=1.5)
    ya = C.verify_axis("log10", y_marks, max_residual_px=1.5)
    return xa, ya, int(ht[2]) + 4, int(hb[1]) - 4, int(vl0[2]) + 4, int(vr0[1]) - 4


def _largest_span_component(mask: np.ndarray, min_span: int = 400) -> np.ndarray:
    lab, n = ndimage.label(mask)
    keep = np.zeros_like(mask)
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        if len(cc) > 0 and (cc.max() - cc.min()) >= min_span:
            keep[lab == i] = True
    return keep


def read_series(image: np.ndarray, at: list[float], which: str) -> list[dict]:
    xa, ya, r0, r1, c0, c1 = calibrate(image)
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)

    if which == "small_w":
        # Red solid curve (1 + omega)
        raw = (r - g > 100) & (r - b > 100)
        mask = np.zeros_like(raw)
        mask[r0:r1, c0:c1] = raw[r0:r1, c0:c1]
        mask = _largest_span_component(mask, min_span=500)
        return [
            C.sample_solid_curve_y(mask, xa, ya, float(x), r0, r1, note_prefix="FigS2D 1+w: ")
            for x in at
        ]

    if which == "large_w":
        # Green solid curve (2 * omega^(1/2)); exclude green tactoid icon and "2 w^(1/2)" label
        raw = (g - r > 60) & (g - b > 60)
        mask = np.zeros_like(raw)
        mask[r0:r1, c0:c1] = raw[r0:r1, c0:c1]
        mask = _largest_span_component(mask, min_span=400)
        return [
            C.sample_solid_curve_y(mask, xa, ya, float(x), r0, r1, note_prefix="FigS2D 2w^(1/2): ")
            for x in at
        ]

    if which == "numerical_min":
        # Black filled circles from numerical minimization (isolated for omega <= 1.15, merged chain above)
        raw = (r < 60) & (g < 60) & (b < 60)
        # Exclude inward ticks along all 4 spines (within 38 px of spines)
        mask = np.zeros_like(raw)
        mask[r0 + 38 : r1 - 38, c0 + 38 : c1 - 38] = raw[r0 + 38 : r1 - 38, c0 + 38 : c1 - 38]
        # Keep only circular markers (area 400..1400 px) or the merged high-omega chain (> 1400 px)
        lab, _ = ndimage.label(mask)
        dots = []
        chain = np.zeros_like(mask)
        for b in C.find_blobs(mask, min_area=400):
            if b["area"] <= 1400:
                dots.append((b["c"], b["r"], float(b["h"])))
            else:
                sl = (slice(b["r0"], b["r1"] + 1), slice(b["c0"], b["c1"] + 1))
                chain[sl] |= lab[sl] == b["label"]

        # Measure true dot radius from isolated circular dots (unoccluded by green curve)
        med_diam = float(np.median([d[2] for d in dots])) if dots else 30.0
        dot_rad = 0.5 * med_diam

        out = []
        for x in at:
            xf = float(x)
            cf = float(xa.to_pixel(xf))
            near_dots = [d for d in dots if abs(d[0] - cf) <= 22.0]
            if near_dots:
                best = min(near_dots, key=lambda d: abs(d[0] - cf))
                c_dot, r_dot, diam = best
                x_dot = float(xa.to_value(c_dot))
                # Estimate local pixel slope from neighboring dots to propagate horizontal dot uncertainty
                nbrs = sorted(dots, key=lambda d: abs(d[0] - c_dot))[:3]
                slope_dot = 0.0
                if len(nbrs) >= 2:
                    slope_dot = float(np.polyfit([d[0] for d in nbrs], [d[1] for d in nbrs], 1)[0])
                half_px = 3.0 + ya.residual_px + abs(slope_dot) * (2.5 + xa.residual_px)
                yv = float(ya.to_value(r_dot))
                ya_, yb_ = float(ya.to_value(r_dot - half_px)), float(ya.to_value(r_dot + half_px))
                out.append(
                    {
                        "x": x_dot,
                        "y": yv,
                        "lo": min(ya_, yb_),
                        "hi": max(ya_, yb_),
                        "px": [c_dot, r_dot],
                        "px_lo": float(r_dot - half_px),
                        "px_hi": float(r_dot + half_px),
                        "note": f"FigS2D dot centroid (d={diam:.0f}px), cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
                    }
                )
            else:
                # In the merged chain (omega >= 1.2), the green 2*omega^(1/2) line runs along the
                # lower edge of the black dots, slightly occluding their bottom pixels. Reconstruct
                # the circle centerline using the unobstructed top edge + dot radius, averaging
                # across a small horizontal window of columns to smooth dot-to-dot scalloping.
                c_int = int(round(cf))
                top_rows = []
                for dc in range(-14, 15):
                    cc = c_int + dc
                    if 0 <= cc < chain.shape[1]:
                        rs = C.runs_1d(chain[r0 + 38 : r1 - 38, cc])
                        if rs:
                            seg = max(rs, key=lambda t: t[1] - t[0])
                            top_rows.append((float(cc) + 0.5, float(seg[0] + r0 + 38)))
                if len(top_rows) < 3:
                    out.append({"x": xf, "y": None, "found": False, "note": f"no chain at col {c_int}"})
                    continue
                xs_t = np.array([p[0] for p in top_rows]) - cf
                ys_t = np.array([p[1] for p in top_rows])
                # Linear fit to local upper envelope of dot chain; add dot_rad * sqrt(1 + slope^2)
                slope, r_top_at_cf = np.polyfit(xs_t, ys_t, 1)
                # Offset vertically from tangent of circle envelope to center: R * sqrt(1 + slope^2)
                # Note upper envelope oscillates slightly between circle crowns (at R*sqrt(1+m^2)) and cusps;
                # use the 25th percentile residual (near crowns) for accurate crown tangent!
                crown_shift = float(np.percentile(ys_t - (slope * xs_t + r_top_at_cf), 20))
                r_center = float(r_top_at_cf + crown_shift + dot_rad * np.hypot(1.0, slope))
                half_px = 0.35 * dot_rad + ya.residual_px + abs(slope) * (0.5 + xa.residual_px)
                yv = float(ya.to_value(r_center))
                ya_, yb_ = float(ya.to_value(r_center - half_px)), float(ya.to_value(r_center + half_px))
                out.append(
                    {
                        "x": xf,
                        "y": yv,
                        "lo": min(ya_, yb_),
                        "hi": max(ya_, yb_),
                        "px": [cf, r_center],
                        "px_lo": float(r_center - half_px),
                        "px_hi": float(r_center + half_px),
                        "note": f"FigS2D merged chain upper-crown+R (R={dot_rad:.1f}px), cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
                    }
                )
        return out

    raise ValueError(which)
