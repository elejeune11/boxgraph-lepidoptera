"""Shared calibrated reader for Main Fig. 6f (4 colored dot series: Droplet 1 & 2, Actin & VASP)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage
import _calib as C

SERIES_NAMES = ("d1_actin", "d1_vasp", "d2_actin", "d2_vasp")


def read_series(image: np.ndarray, at: list[float], which: str) -> list[dict]:
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    sat = np.maximum(np.maximum(r, g), b) - np.minimum(np.minimum(r, g), b)
    dark, gray, hl, vl = C.find_l_spines(image, row_frac=0.4, col_frac=0.4)
    hb, vl0 = hl[-1], vl[0]
    dark_blk = dark & (sat < 30.0)
    xt = C.match_ticks_regular(
        C.find_outward_ticks(dark_blk, hb, (int(vl0[1]) - 10, image.shape[1]), "x", side=+1, band=(6, 25)),
        6,
        tol=5.0,
    )
    # Use outer tick band (22..42 px left of spine) so the green dot at t=0, y=1.0 never occludes tick 1.0
    yt = C.match_ticks_regular(
        C.find_outward_ticks(dark_blk, vl0, (0, int(hb[2]) + 10), "y", side=-1, band=(22, 42), max_width=20),
        8,
        tol=6.0,
    )
    xa = C.verify_axis("linear", list(zip(xt, [0.0, 50.0, 100.0, 150.0, 200.0, 250.0])), max_residual_px=3.0)
    ya = C.verify_axis("linear", list(zip(yt, [3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.0])), max_residual_px=4.0)

    # Sample the 4 dot RGB colors directly from the legend on the right (x in [215, 238])
    c0, c1 = int(xa.to_pixel(215.0)), int(xa.to_pixel(238.0))
    r0, r1 = int(ya.to_pixel(2.7)), int(ya.to_pixel(0.4))
    bands = C.clusters(np.where((sat[r0:r1, c0:c1] > 35.0).any(axis=1))[0] + r0, gap=15)
    if len(bands) != 4:
        raise RuntimeError(f"Expected 4 legend dots in Fig. 6f, found {len(bands)}")
    colors = {}
    for name, (_, b0, b1) in zip(SERIES_NAMES, bands):
        bp = image[b0 : b1 + 1, c0:c1].astype(float)
        bs = sat[b0 : b1 + 1, c0:c1]
        colors[name] = np.median(bp[bs >= 0.85 * float(bs.max())], axis=0)

    rgb = image.astype(float)
    dists = {name: np.linalg.norm(rgb - col[None, None, :], axis=2) for name, col in colors.items()}
    dt = dists[which]
    do = np.minimum.reduce([d for name, d in dists.items() if name != which])
    mask = (sat > 30.0) & (dt < 45.0) & (dt + 15.0 < do)
    # Restrict to plot data area (x <= 205 s, excluding legend on right)
    mask[:, int(xa.to_pixel(205.0)) :] = False
    mask[: int(yt[0]) - 20, :] = False
    mask[int(hb[1]) - 3 :, :] = False

    # Extract connected marker blobs (each dot is ~36x36 px, area ~800-1200 px) using fast find_objects
    lab, n = ndimage.label(mask)
    slices = ndimage.find_objects(lab)
    dots = []
    for i, sl in enumerate(slices, start=1):
        if sl is None:
            continue
        sub = lab[sl] == i
        if int(sub.sum()) < 120:
            continue
        rr_l, cc_l = np.nonzero(sub)
        rr = rr_l + sl[0].start
        cc = cc_l + sl[1].start
        w_b = int(cc.max() - cc.min() + 1)
        if w_b > 55:
            n_sub = max(2, int(round(w_b / 38.0)))
            edges = np.linspace(cc.min(), cc.max() + 1, n_sub + 1)
            for s in range(n_sub):
                sel = (cc >= edges[s]) & (cc < edges[s + 1])
                if sel.sum() >= 80:
                    rc, cc_m = float(rr[sel].mean()) + 0.5, float(cc[sel].mean()) + 0.5
                    dots.append((float(xa.to_value(cc_m)), float(ya.to_value(rc)), cc_m, rc))
        else:
            rc, cc_m = float(rr.mean()) + 0.5, float(cc.mean()) + 0.5
            dots.append((float(xa.to_value(cc_m)), float(ya.to_value(rc)), cc_m, rc))
    dots.sort(key=lambda d: d[0])

    h_px = 4.0 + ya.residual_px
    out = []
    if at:
        for tx in at:
            near = min(dots, key=lambda d: abs(d[0] - tx))
            if abs(near[0] - tx) > 6.0:
                out.append({"x": float(tx), "y": None, "found": False, "note": f"no {which} dot near t={tx} s"})
                continue
            xv, yv, cc_m, rc = near
            ya_, yb_ = float(ya.to_value(rc - h_px)), float(ya.to_value(rc + h_px))
            out.append(
                {
                    "x": float(round(xv, 2)),
                    "y": yv,
                    "lo": min(ya_, yb_),
                    "hi": max(ya_, yb_),
                    "px": [cc_m, rc],
                    "px_lo": float(rc - h_px),
                    "px_hi": float(rc + h_px),
                    "note": f"Fig6f {which} dot at t~{xv:.1f}s, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
                }
            )
    else:
        for xv, yv, cc_m, rc in dots:
            ya_, yb_ = float(ya.to_value(rc - h_px)), float(ya.to_value(rc + h_px))
            out.append(
                {
                    "x": float(round(xv, 2)),
                    "y": yv,
                    "lo": min(ya_, yb_),
                    "hi": max(ya_, yb_),
                    "px": [cc_m, rc],
                    "px_lo": float(rc - h_px),
                    "px_hi": float(rc + h_px),
                    "note": f"Fig6f {which}: cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
                }
            )
    return out
