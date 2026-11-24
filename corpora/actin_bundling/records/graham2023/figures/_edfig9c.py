"""Shared calibrated reader for Extended Data Fig. 9c droplet size distributions (3%, 5%, 10%, 15% PEG)."""

from __future__ import annotations

import numpy as np
import _calib as C

PEG_NAMES = ("3peg", "5peg", "10peg", "15peg")


def read_dist(image: np.ndarray, at: list[float], which: str) -> list[dict]:
    dark, gray, hl, vl = C.find_l_spines(image, row_frac=0.25, col_frac=0.4)
    vl0 = vl[0]
    yt = C.find_outward_ticks(dark, vl0, (0, image.shape[0]), "y", side=-1, band=(6, 30), max_width=24)
    yt = C.match_ticks_regular(yt, 8, tol=6.0)
    ya = C.verify_axis("linear", list(zip(yt, [35.0, 30.0, 25.0, 20.0, 15.0, 10.0, 5.0, 0.0])), max_residual_px=5.0)
    r_zero = int(round(yt[-1]))
    xt = C.find_outward_ticks(dark, (float(r_zero), r_zero - 4, r_zero + 4), (int(vl0[1]) - 15, image.shape[1]), "x", side=+1, band=(10, 38), max_width=24)
    xt = C.match_ticks_regular(xt, 5, tol=6.0)
    xa = C.verify_axis("linear", list(zip(xt, [0.0, 1.0, 2.0, 3.0, 4.0])), max_residual_px=5.0)

    # Sample the 4 curve RGBs directly from the upper-right legend
    c0, c1 = int(xa.to_pixel(2.2)), int(xa.to_pixel(2.75))
    r0, r1 = int(ya.to_pixel(35.0)), int(ya.to_pixel(18.0))
    patch = image[r0:r1, c0:c1].astype(float)
    sat = np.maximum(np.maximum(patch[:, :, 0], patch[:, :, 1]), patch[:, :, 2]) - np.minimum(
        np.minimum(patch[:, :, 0], patch[:, :, 1]), patch[:, :, 2]
    )
    bands = C.clusters(np.where((sat > 50).mean(axis=1) > 0.4)[0] + r0, gap=10)
    if len(bands) != 4:
        raise RuntimeError(f"Expected 4 legend bands in ED Fig. 9c, found {len(bands)}")
    colors = {}
    for name, (_, b0, b1) in zip(PEG_NAMES, bands):
        bp = image[b0 : b1 + 1, c0 + 5 : c1 - 5].astype(float)
        bs = np.maximum(np.maximum(bp[:, :, 0], bp[:, :, 1]), bp[:, :, 2]) - np.minimum(
            np.minimum(bp[:, :, 0], bp[:, :, 1]), bp[:, :, 2]
        )
        colors[name] = np.median(bp[bs >= 0.8 * float(bs.max())], axis=0)

    rgb = image.astype(float)
    sat_all = np.maximum(np.maximum(rgb[:, :, 0], rgb[:, :, 1]), rgb[:, :, 2]) - np.minimum(
        np.minimum(rgb[:, :, 0], rgb[:, :, 1]), rgb[:, :, 2]
    )
    dists = {name: np.linalg.norm(rgb - col[None, None, :], axis=2) for name, col in colors.items()}
    dt = dists[which]
    do = np.minimum.reduce([d for name, d in dists.items() if name != which])
    mask = (sat_all > 40.0) & (dt < 55.0) & (dt + 15.0 < do)
    # Mask upper-right legend
    mask[: int(ya.to_pixel(16.0)), int(xa.to_pixel(2.1)) :] = False

    row_min = int(yt[0]) - 20
    row_max = r_zero + 15
    out = []
    for x in at:
        res = C.sample_solid_curve_y(mask, xa, ya, x, row_min, row_max, search_dx=4, note_prefix=f"ED9c {which}: ")
        if (not res.get("found", res.get("y") is not None)) and x >= 2.3:
            # Overlapping zero-baseline tail (e.g. yellow 15% PEG covers red/blue at y=0 for D >= 2.4 um)
            cf = float(xa.to_pixel(x))
            c = int(round(cf))
            r_zero = int(round(ya.to_pixel(0.0)))
            band = sat_all[r_zero - 18 : r_zero + 10, max(0, c - 4) : c + 5] > 40.0
            if band.any():
                half_px = 6.0 + ya.residual_px
                yv = 0.0
                ya_, yb_ = float(ya.to_value(r_zero - half_px)), float(ya.to_value(r_zero + half_px))
                res = {
                    "x": float(x),
                    "y": yv,
                    "lo": max(0.0, min(ya_, yb_)),
                    "hi": max(ya_, yb_),
                    "px": [cf, float(r_zero)],
                    "px_lo": float(r_zero - half_px),
                    "px_hi": float(r_zero + half_px),
                    "note": f"ED9c {which}: overlapping zero baseline, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
                }
        out.append(res)
    return out
