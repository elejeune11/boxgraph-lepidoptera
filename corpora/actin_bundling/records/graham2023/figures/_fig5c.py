"""Shared calibrated reader for Main Fig. 5c (Bending and Surface energy vs Aspect ratio)."""

from __future__ import annotations

import numpy as np
import _calib as C


def calibrate(image: np.ndarray):
    dark, gray, hl, vl = C.find_l_spines(image, row_frac=0.5, col_frac=0.5)
    hb = hl[-1]
    vl0 = vl[0]
    xt = C.find_outward_ticks(dark, hb, (int(vl0[2]), image.shape[1]), "x", side=+1, band=(6, 25))
    yt = C.find_outward_ticks(dark, vl0, (0, int(hb[1])), "y", side=-1, band=(6, 25))
    xt = C.match_ticks_regular(xt, 6, tol=4.0)
    yt = C.match_ticks_regular(yt, 3, tol=3.0)
    xa = C.verify_axis("linear", list(zip(xt, [1.0, 3.0, 5.0, 7.0, 9.0, 11.0])), max_residual_px=2.0)
    ya = C.verify_axis("linear", list(zip(yt, [20.0, 10.0, 0.0])), max_residual_px=1.0)
    row_min = int(yt[0]) - 25
    row_max = int(hb[1]) - 5
    return xa, ya, row_min, row_max, int(vl0[2])


def read_curve(image: np.ndarray, at: list[float], which: str) -> list[dict]:
    xa, ya, row_min, row_max, c_left = calibrate(image)
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)

    # Dark green Bending curve: RGB ~(24, 142, 71)
    bend = (g - r > 45) & (g > 70) & (r < 110)
    # Mask legend box in upper-left of plot area
    bend[: int(ya.to_pixel(14.0)), : int(xa.to_pixel(7.2))] = False

    if which == "bending":
        mask = bend
        prefix = "Fig5c bending: "
    elif which == "surface":
        # Pale grey-green Surface curve: RGB ~(176, 209, 207)
        dist = np.sqrt((r - 176.0) ** 2 + (g - 209.0) ** 2 + (b - 207.0) ** 2)
        mask = (dist < 42.0) & ~bend
        mask[: int(ya.to_pixel(14.0)), : int(xa.to_pixel(7.2))] = False
        prefix = "Fig5c surface: "
    else:
        raise ValueError(which)

    out = []
    for x in at:
        if which == "bending" and 9.45 <= x <= 9.65:
            # Near-vertical drop: take full vertical stroke span in column directly as honest interval
            c_f = float(xa.to_pixel(x))
            c = int(round(c_f))
            rs = [(a + row_min, b_ + row_min) for a, b_ in C.runs_1d(mask[row_min : row_max + 1, c], gap=3)]
            if not rs:
                out.append({"x": float(x), "y": None, "found": False, "note": f"no bend pixel at col {c}"})
                continue
            r0 = min(a for a, _ in rs)
            r1 = max(b_ for _, b_ in rs)
            r_mid = 0.5 * (r0 + r1) + 0.5
            half_px = 0.5 * (r1 - r0 + 1) + 1.0 + ya.residual_px
            yv = float(ya.to_value(r_mid))
            ya_, yb_ = float(ya.to_value(r_mid - half_px)), float(ya.to_value(r_mid + half_px))
            out.append(
                {
                    "x": float(x),
                    "y": yv,
                    "lo": min(ya_, yb_),
                    "hi": max(ya_, yb_),
                    "px": [c_f, float(r_mid)],
                    "px_lo": float(r_mid - half_px),
                    "px_hi": float(r_mid + half_px),
                    "note": f"{prefix}steep drop span={r1-r0+1}px, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
                }
            )
        else:
            out.append(C.sample_solid_curve_y(mask, xa, ya, x, row_min, row_max, note_prefix=prefix))
    return out
