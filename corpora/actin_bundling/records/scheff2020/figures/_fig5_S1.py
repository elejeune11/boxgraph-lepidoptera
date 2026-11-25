"""Shared readers for Fig. 5(c), Fig. 5(d), and ESI Fig. S1 (left & right)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage
import _helpers as H


def _cal_fig5c(img: np.ndarray):
    box = H.find_box(img)
    xt = H.find_ticks(img, box, "top", min_len=10)
    if len(xt) != 7:
        raise ValueError(f"fig5c: expected 7 top x ticks, got {len(xt)}: {xt}")
    xpairs = [(box["c0"], 0.0)] + [(xt[i][0], 5.0 * (i + 1)) for i in range(7)] + [(box["c1"], 40.0)]
    yt = H.find_ticks(img, box, "right", min_len=10)
    if len(yt) != 4:
        raise ValueError(f"fig5c: expected 4 right y ticks, got {len(yt)}: {yt}")
    ypairs = [(yt[i][0], 2.0 - 0.25 * i) for i in range(4)] + [(box["r1"], 1.0)]
    return box, H.AxisCal(xpairs, max_res_px=3.0), H.AxisCal(ypairs, max_res_px=3.0)


def _cal_figS1_left(img: np.ndarray):
    box = H.find_box(img)
    xt = H.find_ticks(img, box, "top", min_len=10)
    if len(xt) != 3:
        raise ValueError(f"figS1_left: expected 3 top x ticks, got {len(xt)}: {xt}")
    xpairs = [(box["c0"], 0.0), (xt[0][0], 10.0), (xt[1][0], 20.0), (xt[2][0], 30.0), (box["c1"], 40.0)]
    yt = H.find_ticks(img, box, "right", min_len=10)
    if len(yt) != 5:
        raise ValueError(f"figS1_left: expected 5 right y ticks, got {len(yt)}: {yt}")
    ypairs = [(yt[i][0], 1.4 - 0.1 * i) for i in range(5)]
    return box, H.AxisCal(xpairs, max_res_px=2.5), H.AxisCal(ypairs, max_res_px=2.5)


def _cal_figS1_right(img: np.ndarray):
    box = H.find_box(img)
    xt = H.find_ticks(img, box, "top", min_len=10)
    if len(xt) != 5:
        raise ValueError(f"figS1_right: expected 5 top x ticks, got {len(xt)}: {xt}")
    xpairs = [(box["c0"], 0.0)] + [(xt[i][0], 10.0 * (i + 1)) for i in range(5)]
    yt = H.find_ticks(img, box, "right", min_len=10)
    if len(yt) != 8:
        raise ValueError(f"figS1_right: expected 8 right y ticks, got {len(yt)}: {yt}")
    ypairs = [(yt[i][0], 4.5 - 0.5 * i) for i in range(8)]
    return box, H.AxisCal(xpairs, max_res_px=2.5), H.AxisCal(ypairs, max_res_px=2.5)


_CALS = {
    "fig5c": _cal_fig5c,
    "figS1_left": _cal_figS1_left,
    "figS1_right": _cal_figS1_right,
}


def read_tactoid_curve(img: np.ndarray, at: list[float], panel: str, which: str) -> list[dict]:
    """Read one of the three theoretical curves ('best', 'min', 'max') in Fig. 5c or Fig. S1."""
    box, xa, ya = _CALS[panel](img)
    r = img[:, :, 0].astype(np.int16)
    g = img[:, :, 1].astype(np.int16)
    b = img[:, :, 2].astype(np.int16)
    if which == "best":
        mask = (r > 150) & (r - g > 80) & (r - b > 80)
    elif which == "min":
        mask = (b > 140) & (g > 100) & (b - r > 70) & (g - r > 50)
    elif which == "max":
        mask = (b > 65) & (b - g > 20) & (b - r > 12) & (g < 95) & (r < 115)
    else:
        raise ValueError(f"Unknown curve {which!r}")
    leg = H.find_inner_frame_box(img, box)
    if leg is not None:
        c0, r0, c1, r1 = leg
        mask[r0:r1 + 1, c0:c1 + 1] = False
    return H.sample_curve_at(mask, xa, ya, at, box, fit_half_win=65, max_gap_px=140)


def read_tactoid_data(img: np.ndarray, panel: str) -> list[dict]:
    """Read black square markers and their vertical ±1 SD error bars in Fig. 5c or Fig. S1."""
    box, xa, ya = _CALS[panel](img)
    r0_off, c0_off = box["r0_e"] + 12, box["c0_e"] + 6
    r1_lim, c1_lim = box["r1_s"] - 6, box["c1_s"] - 6
    sub = img[r0_off:r1_lim, c0_off:c1_lim]

    mx = sub.max(axis=2).astype(np.int16)
    mn = sub.min(axis=2).astype(np.int16)
    # Combined ink catches squares even where colored curves cross over them
    ink = mn < 135
    # Pure neutral dark ink distinguishes squares and error-bar stems from colored curves
    black = (mx < 105) & ((mx - mn) < 35)
    # Stem mask: neutral dark/grey pixels for tracing vertical error bars
    full_mx = img.max(axis=2).astype(np.int16)
    full_mn = img.min(axis=2).astype(np.int16)
    stem_mask = (full_mx < 175) & ((full_mx - full_mn) < 38)
    # Exclude top/bottom/left/right spines from stem mask
    stem_mask[:box["r0_e"] + 3, :] = False
    stem_mask[box["r1_s"] - 2:, :] = False
    stem_mask[:, :box["c0_e"] + 3] = False
    stem_mask[:, box["c1_s"] - 2:] = False

    leg = H.find_inner_frame_box(img, box)
    if leg is not None:
        lc0, lr0, lc1, lr1 = leg
        ink[max(0, lr0 - r0_off):max(0, lr1 - r0_off + 1), max(0, lc0 - c0_off):max(0, lc1 - c0_off + 1)] = False
        stem_mask[lr0:lr1 + 1, lc0:lc1 + 1] = False

    op = H.open_rect(ink, 22, 22)
    out = []
    for _, rsl, csl, loc in H.labeled_components(op):
        h = rsl.stop - rsl.start
        w = csl.stop - csl.start
        if w < 26 or w > 58 or h < 26 or h > 72:
            continue
        bf = float(black[rsl, csl][loc].mean())
        if bf < 0.35:
            continue
        # Refine center using neutral black pixels inside the bounding box so a curve grazing
        # one edge of the square does not bias the marker centroid
        loc_blk = black[rsl, csl] & loc
        if loc_blk.sum() >= 50:
            rr, cc = np.nonzero(loc_blk)
            r_center = r0_off + rsl.start + float(rr.mean()) + 0.5
            c_center = c0_off + csl.start + float(cc.mean()) + 0.5
            r_top = r0_off + rsl.start + int(rr.min())
            r_bot = r0_off + rsl.start + int(rr.max())
        else:
            rr, cc = np.nonzero(loc)
            r_center = r0_off + rsl.start + float(rr.mean()) + 0.5
            c_center = c0_off + csl.start + float(cc.mean()) + 0.5
            r_top = r0_off + rsl.start
            r_bot = r0_off + rsl.stop - 1

        # Trace vertical error bar upward and downward from the square
        # Allow gap up to 24 px across colored curves crossing the error bar stem
        r_up, r_dn = H.trace_vertical_ebar(
            stem_mask, c_center, r_top - 2, r_bot + 2,
            r_min_bound=box["r0_e"] + 3, r_max_bound=box["r1_s"] - 3,
            col_half_win=3, max_gap=24,
        )
        # Require stem to extend at least 10 px beyond marker edge to count as a real error bar
        if r_up is not None and (r_top - r_up) < 10.0:
            r_up = None
        if r_dn is not None and (r_dn - r_bot) < 10.0:
            r_dn = None
        # Clip check: if one bar end is very close to bottom spine and upper bar exists, keep bottom spine as lower bound
        has_bar = (r_up is not None) or (r_dn is not None)
        note = "square with ±1 SD bar" if has_bar else "single-droplet square (no error bar)"
        if r_up is not None and r_dn is None and (box["r1_s"] - r_bot) < (r_top - r_up):
            # Lower bar extends to/past bottom axis frame
            r_dn = float(box["r1_s"] - 2)
            note += " (lower end clipped at bottom axis)"
        rd = H.make_marker_reading(
            xa, ya, c_center, r_center, r_up, r_dn,
            marker_half_h=0.5 * (r_bot - r_top + 1),
            note=note, include_ebar=has_bar,
        )
        out.append(rd)
    return sorted(out, key=lambda d: d["x"])


def read_fig5d(img: np.ndarray, series: str) -> list[dict]:
    """Read Fig. 5(d): series='K_over_gamma' (black squares, left axis) or 'omega' (red circles, right axis)."""
    box = H.find_box(img)
    xt = H.find_ticks(img, box, "top", min_len=10)
    if len(xt) != 4:
        raise ValueError(f"fig5d: expected 4 top x ticks, got {len(xt)}: {xt}")
    xpairs = [(xt[0][0], 100.0), (xt[1][0], 150.0), (xt[2][0], 200.0), (xt[3][0], 250.0)]
    xa = H.AxisCal(xpairs, max_res_px=2.5)

    r_ch = img[:, :, 0].astype(np.int16)
    g_ch = img[:, :, 1].astype(np.int16)
    b_ch = img[:, :, 2].astype(np.int16)
    mx = np.maximum(np.maximum(r_ch, g_ch), b_ch)
    mn = np.minimum(np.minimum(r_ch, g_ch), b_ch)

    if series == "K_over_gamma":
        lt = H.find_ticks(img, box, "left", min_len=10)
        if len(lt) != 6:
            raise ValueError(f"fig5d: expected 6 left y ticks, got {len(lt)}: {lt}")
        ypairs = [(box["r0"], 2.5)] + [(lt[i][0], 2.25 - 0.25 * i) for i in range(6)]
        ya = H.AxisCal(ypairs, max_res_px=2.5)
        raw_mask = (mx < 100) & ((mx - mn) < 35)
    elif series == "omega":
        rt = H.find_ticks(img, box, "right", min_len=10)
        if len(rt) != 6:
            raise ValueError(f"fig5d: expected 6 right y ticks, got {len(rt)}: {rt}")
        ypairs = [(rt[i][0], 3.0 - 0.5 * i) for i in range(6)] + [(box["r1"], 0.0)]
        ya = H.AxisCal(ypairs, max_res_px=2.5)
        raw_mask = (r_ch > 150) & (r_ch - g_ch > 80) & (r_ch - b_ch > 80)
    else:
        raise ValueError(f"Unknown fig5d series {series!r}")

    # Restrict strictly inside plot box and exclude legend
    mask = np.zeros_like(raw_mask)
    r0_i, r1_i = box["r0_e"] + 10, box["r1_s"] - 10
    c0_i, c1_i = box["c0_e"] + 10, box["c1_s"] - 10
    mask[r0_i:r1_i, c0_i:c1_i] = raw_mask[r0_i:r1_i, c0_i:c1_i]
    leg = H.find_inner_frame_box(img, box)
    if leg is not None:
        lc0, lr0, lc1, lr1 = leg
        mask[lr0:lr1 + 1, lc0:lc1 + 1] = False

    # Seal narrow vertical crossings (black stem crossing red circle rim at L=90 nm) without merging horizontal caps
    sealed = H.close_rect(mask, 3, 15) if series == "omega" else mask
    filled = ndimage.binary_fill_holes(sealed)
    markers_only = H.open_rect(filled, 30, 30)
    out = []
    for _, rsl, csl, loc in H.labeled_components(markers_only):
        h = rsl.stop - rsl.start
        w = csl.stop - csl.start
        if not (35 <= w <= 95 and 35 <= h <= 95):
            continue
        rr, cc = np.nonzero(loc)
        r_center = rsl.start + float(rr.mean()) + 0.5
        c_center = csl.start + float(cc.mean()) + 0.5
        r_top = rsl.start
        r_bot = rsl.stop - 1
        # Trace range error bar upward and downward in same color mask; bridge black cap crossings (up to 32 px)
        r_up, r_dn = H.trace_vertical_ebar(
            mask, c_center, r_top - 1, r_bot + 1,
            r_min_bound=r0_i, r_max_bound=r1_i,
            col_half_win=8, max_gap=32,
        )
        # Keep even short range bars (e.g. omega=0.20 -> 0.15 or 0.45 -> 0.38) if >= 2 px beyond rim
        min_ext = 2.0 if series == "omega" else 6.0
        if r_up is not None and (r_top - r_up) < min_ext:
            r_up = None
        if r_dn is not None and (r_dn - r_bot) < min_ext:
            r_dn = None
        has_bar = (r_up is not None) or (r_dn is not None)
        note = f"{series} with min/max bounding fit range bar"
        if series == "omega" and r_up is not None and r_dn is None:
            r_dn = float(r_bot)
            note += " (lower end hidden inside marker)"
        rd = H.make_marker_reading(
            xa, ya, c_center, r_center, r_up, r_dn,
            marker_half_h=0.5 * h,
            note=note,
            include_ebar=has_bar,
        )
        out.append(rd)
    return sorted(out, key=lambda d: d["x"])
