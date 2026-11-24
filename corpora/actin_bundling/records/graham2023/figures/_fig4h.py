"""Shared calibrated reader for Main Fig. 4h (predicted critical T_ring vs D for 7 gamma_I values,
plus the overlaid black solid power-law fit to experimental AR < 1.1 droplets).
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage
import _calib as C

# Legend order top-to-bottom in Fig. 4h: 15.0, 12.5, 10.0, 7.5, 5.0, 2.5, 1.0 pN
GI_NAMES = ("gI15", "gI12p5", "gI10", "gI7p5", "gI5", "gI2p5", "gI1")


def calibrate(image: np.ndarray):
    dark, gray, hl, vl = C.find_l_spines(image, row_frac=0.4, col_frac=0.4)
    ht, hb = hl[0], hl[-1]
    vl0, vr0 = vl[0], vl[-1]
    xt = C.find_outward_ticks(dark, hb, (int(vl0[1]) - 10, int(vr0[2]) + 10), "x", side=+1, band=(6, 25))
    yt = C.find_outward_ticks(dark, vl0, (int(ht[1]) - 10, int(hb[2]) + 10), "y", side=-1, band=(6, 25))
    if len(xt) != 5:
        raise RuntimeError(f"Expected 5 x-ticks (0.5,1,2,3,4) in Fig. 4h, found {xt}")
    yt = C.match_ticks_regular(yt, 5, tol=4.0)
    xa = C.verify_axis("linear", list(zip(xt, [0.5, 1.0, 2.0, 3.0, 4.0])), max_residual_px=2.5)
    ya = C.verify_axis("linear", list(zip(yt, [2.0, 1.5, 1.0, 0.5, 0.0])), max_residual_px=2.5)

    # Sample the 7 gamma_I RGB colors directly from the legend lines inside the plot
    r0, r1 = int(ya.to_pixel(1.85)), int(ya.to_pixel(0.85))
    c0, c1 = int(xa.to_pixel(0.65)), int(xa.to_pixel(1.2))
    patch = image[r0:r1, c0:c1].astype(float)
    sat = np.maximum(np.maximum(patch[:, :, 0], patch[:, :, 1]), patch[:, :, 2]) - np.minimum(
        np.minimum(patch[:, :, 0], patch[:, :, 1]), patch[:, :, 2]
    )
    bands = C.clusters(np.where((sat > 40).mean(axis=1) > 0.4)[0] + r0, gap=10)
    if len(bands) != 7:
        raise RuntimeError(f"Expected 7 legend bands in Fig. 4h, found {len(bands)}")
    colors = {}
    for name, (_, b0, b1) in zip(GI_NAMES, bands):
        bp = image[b0 : b1 + 1, c0 + 10 : c1 - 10].astype(float)
        bs = np.maximum(np.maximum(bp[:, :, 0], bp[:, :, 1]), bp[:, :, 2]) - np.minimum(
            np.minimum(bp[:, :, 0], bp[:, :, 1]), bp[:, :, 2]
        )
        colors[name] = np.median(bp[bs >= 0.8 * float(bs.max())], axis=0)

    bounds = (int(ht[2]) + 3, int(hb[1]) - 3, int(vl0[2]) + 3, int(vr0[1]) - 3)
    leg_box = (int(ht[2]) + 3, int(ya.to_pixel(0.85)), int(vl0[2]) + 3, int(xa.to_pixel(1.85)))
    return xa, ya, colors, bounds, leg_box


def _norm_chroma(arr: np.ndarray) -> np.ndarray:
    """White- and black-blend invariant unit chroma direction (rgb - min) / ||rgb - min||."""
    c = arr - np.min(arr, axis=-1, keepdims=True)
    n = np.linalg.norm(c, axis=-1, keepdims=True)
    return c / np.maximum(n, 1e-6)


def read_gi(image: np.ndarray, at: list[float], which: str) -> list[dict]:
    """Read circular markers and vertical s.d. error bars for one gamma_I series in Fig. 4h."""
    xa, ya, colors, (rmin, rmax, cmin, cmax), (lr0, lr1, lc0, lc1) = calibrate(image)
    rgb = image.astype(float)
    sat = np.maximum(np.maximum(rgb[:, :, 0], rgb[:, :, 1]), rgb[:, :, 2]) - np.minimum(
        np.minimum(rgb[:, :, 0], rgb[:, :, 1]), rgb[:, :, 2]
    )
    u_img = _norm_chroma(rgb)
    u_cols = {name: _norm_chroma(col) for name, col in colors.items()}
    dists = {name: np.linalg.norm(u_img - u_col[None, None, :], axis=2) for name, u_col in u_cols.items()}
    d_t = dists[which]
    d_o = np.minimum.reduce([d for name, d in dists.items() if name != which])
    # Match pure colored pixels by blend-invariant unit chroma direction
    mask = (sat > 40.0) & (d_t < 0.14) & (d_t + 0.02 < d_o)
    mask[lr0:lr1, lc0:lc1] = False
    mask[:rmin, :] = False
    mask[rmax + 1 :, :] = False
    mask[:, :cmin] = False
    mask[:, cmax + 1 :] = False

    H, W = mask.shape
    out = []
    for x in at:
        cf = float(xa.to_pixel(x))
        c = int(round(cf))
        # Narrow strip +/- 4 px for vertical error bar + marker stem; wider strip +/- 18 px for marker disk
        c_s0, c_s1 = max(cmin, c - 4), min(cmax, c + 4)
        c_m0, c_m1 = max(cmin, c - 18), min(cmax, c + 18)

        # Vertical extent at marker center column
        col_hit = mask[rmin : rmax + 1, c_s0 : c_s1 + 1].any(axis=1)
        rs = [(a + rmin, b + rmin) for a, b in C.runs_1d(col_hit, gap=8)]
        if not rs:
            # Try +/- 10 px horizontally if slightly shifted
            col_hit2 = mask[rmin : rmax + 1, max(cmin, c - 10) : min(cmax, c + 10) + 1].any(axis=1)
            rs = [(a + rmin, b + rmin) for a, b in C.runs_1d(col_hit2, gap=8)]
        if not rs:
            out.append({"x": float(x), "y": None, "found": False, "note": f"no {which} marker near D={x}"})
            continue

        # Pick the vertical run that has the largest total area in the marker strip
        seg = max(rs, key=lambda t: int(mask[t[0] : t[1] + 1, c_m0 : c_m1 + 1].sum()))
        r_top, r_bot = seg

        # Distinguish circular marker (~30-38 px wide) from vertical error bar stem (~10-14 px wide)
        widths = mask[r_top : r_bot + 1, c_m0 : c_m1 + 1].sum(axis=1)
        wmax = int(widths.max())
        # Marker body rows have width near wmax (at least max(16, 0.72*wmax))
        m_rows = np.where(widths >= max(16, int(0.72 * wmax)))[0] + r_top
        if len(m_rows) >= 3:
            r_mid = float(m_rows.mean()) + 0.5
            m_half = 0.5 * float(m_rows.max() - m_rows.min() + 1)
        else:
            r_mid = 0.5 * (r_top + r_bot) + 0.5
            m_half = 17.0

        # Reading uncertainty on the marker center itself
        h_px = 3.0 + ya.residual_px
        yv = float(ya.to_value(r_mid))
        ya_, yb_ = float(ya.to_value(r_mid - h_px)), float(ya.to_value(r_mid + h_px))
        lo, hi = min(ya_, yb_), max(ya_, yb_)

        # Check if vertical error bar protrudes past circular marker (~20 px radius)
        up_ext = r_mid - r_top
        dn_ext = r_bot - r_mid
        has_bar = (up_ext > max(23.0, m_half + 6.0)) or (dn_ext > max(23.0, m_half + 6.0))
        if has_bar:
            # Use measured ends (if one end is clipped/occluded by a neighbor, take symmetric max if very asymmetric)
            r_u = float(r_top)
            r_l = float(r_bot)
            bar_hi = float(ya.to_value(r_u))
            bar_lo = float(ya.to_value(r_l))
            bar_hi_out = float(ya.to_value(r_u - h_px))
            bar_lo_out = float(ya.to_value(r_l + h_px))
            note = f"Fig4h {which}: s.d. bar ({r_top}..{r_bot}), cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px"
            px_lo, px_hi = r_u, r_l
        else:
            # Error bar hidden inside marker: bound by marker radius (~17 px)
            mr = min(max(m_half, 14.0), 20.0)
            bar_hi = float(ya.to_value(r_mid - mr))
            bar_lo = max(0.0, float(ya.to_value(r_mid + mr)))
            bar_hi_out = float(ya.to_value(r_mid - mr - h_px))
            bar_lo_out = max(0.0, float(ya.to_value(r_mid + mr + h_px)))
            note = f"Fig4h {which}: s.d. bar hidden inside marker (r~{mr:.0f}px), cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px"
            px_lo, px_hi = r_mid - mr, r_mid + mr

        out.append(
            {
                "x": float(x),
                "y": yv,
                "lo": lo,
                "hi": hi,
                "bar_lo": min(bar_lo, bar_hi),
                "bar_hi": max(bar_lo, bar_hi),
                "bar_lo_outer": min(bar_lo_out, bar_hi_out),
                "bar_hi_outer": max(bar_lo_out, bar_hi_out),
                "px": [cf, r_mid],
                "px_lo": float(px_lo),
                "px_hi": float(px_hi),
                "note": note,
            }
        )
    return out


def read_exp_fit(image: np.ndarray, at: list[float]) -> list[dict]:
    """Read the overlaid black solid power-law experimental trendline in Fig. 4h."""
    xa, ya, _, (rmin, rmax, cmin, cmax), (lr0, lr1, lc0, lc1) = calibrate(image)
    gray = image.astype(float).mean(axis=2)
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    sat = np.maximum(np.maximum(r, g), b) - np.minimum(np.minimum(r, g), b)

    # Pure black pixels inside the plot box (excluding legend and spines)
    black = (gray < 40.0) & (sat < 20.0)
    black[:rmin, :] = False
    black[rmax + 1 :, :] = False
    black[:, :cmin] = False
    black[:, cmax + 1 :] = False
    black[lr0:lr1, lc0:lc1] = False

    # Scatter squares are ~26x26 solid black; morphological opening with 19x19 isolates squares
    squares = ndimage.binary_opening(black, structure=np.ones((19, 19), dtype=bool))
    squares_dil = ndimage.binary_dilation(squares, iterations=4)
    thin_line = black & ~squares_dil

    # Keep only thin-line segments (vertical run height 5..20 px)
    out = []
    for x in at:
        res = C.sample_solid_curve_y(
            thin_line,
            xa,
            ya,
            x,
            rmin,
            rmax,
            pick_fn=lambda rs, col: min(
                [s for s in rs if 5 <= (s[1] - s[0] + 1) <= 22] or rs,
                key=lambda s: abs((s[1] - s[0] + 1) - 12),
            ),
            search_dx=3,
            note_prefix="Fig4h exp power fit: ",
        )
        out.append(res)
    return out
