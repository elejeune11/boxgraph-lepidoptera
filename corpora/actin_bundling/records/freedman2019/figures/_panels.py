"""Panel-specific figure readers for freedman2019 built on `_calib.py`.

Every function takes `(image, at)` and returns a list of BoxGraph code-mode reading dicts.
"""
from __future__ import annotations

from typing import Sequence

import numpy as np
from scipy import ndimage

import _calib as C


# ---------------------------------------------------------------------------
# Fig. 2C (p. 3): domain length vs Delta l_xl (lattice model & AFINES)
# ---------------------------------------------------------------------------

def read_fig2c(img: np.ndarray, at: Sequence[float], series: str) -> list[dict]:
    # Bottom ticks: 30, 100, 300, 1000 are true log10; leftmost tick (at log10(10)) is labeled "0"
    xa, ya, box, _ = C.calibrate_mathematica_panel(
        img,
        x_tick_vals=[10.0, 30.0, 100.0, 300.0, 1000.0],
        y_tick_vals=[6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0],
        x_scale="log10",
        y_scale="linear",
    )
    c0, r0, c1, r1 = box
    leg_r1 = int(round(ya.to_pixel(4.1)))
    leg_c1 = int(round(xa.to_pixel(150.0)))
    legend_rects = [(r0, leg_r1, c0, leg_c1)]

    x_positions = [
        (0.0, 10.0),
        (27.0, 27.0),
        (50.0, 50.0),
        (100.0, 100.0),
        (300.0, 300.0),
        (500.0, 500.0),
        (1000.0, 1000.0),
    ]
    if series == "lattice":
        mask = C.black_mask(img, thresh=85.0, sat_max=35.0)
        return C.read_mathematica_series_at_x(
            img, xa, ya, box, mask, x_positions,
            legend_rects=legend_rects, cap_half_w=14, filled_rad=14, inset_ticks_px=24,
        )
    elif series == "afines":
        mask = C.red_mask(img, r_min=150.0, diff=50.0)
        return C.read_mathematica_series_at_x(
            img, xa, ya, box, mask, x_positions, legend_rects=legend_rects, cap_half_w=14
        )
    raise ValueError(f"Unknown Fig. 2C series {series!r}")


# ---------------------------------------------------------------------------
# Fig. 2D (p. 3): domain length vs L_p (lattice model & AFINES)
# ---------------------------------------------------------------------------

def read_fig2d(img: np.ndarray, at: Sequence[float], series: str) -> list[dict]:
    xa, ya, box, _ = C.calibrate_mathematica_panel(
        img,
        x_tick_vals=[5.0, 10.0, 50.0, 100.0],
        y_tick_vals=[4.0, 3.5, 3.0, 2.5, 2.0, 1.5],
        x_scale="log10",
        y_scale="linear",
    )
    c0, r0, c1, r1 = box
    # Legend ("lattice model" / "AFINES") sits in y > 3.18, x < 36
    leg_r1 = int(round(ya.to_pixel(3.18)))
    leg_c1 = int(round(xa.to_pixel(36.0)))
    legend_rects = [(r0, leg_r1, c0, leg_c1)]

    x_positions = [(v, v) for v in (2.5, 7.5, 12.5, 17.0, 25.0, 75.0, 125.0)]
    if series == "lattice":
        mask = C.black_mask(img, thresh=85.0, sat_max=35.0)
        return C.read_mathematica_series_at_x(
            img, xa, ya, box, mask, x_positions,
            legend_rects=legend_rects, cap_half_w=14, filled_rad=14, inset_ticks_px=24,
        )
    elif series == "afines":
        mask = C.red_mask(img, r_min=150.0, diff=50.0)
        return C.read_mathematica_series_at_x(
            img, xa, ya, box, mask, x_positions, legend_rects=legend_rects, cap_half_w=14
        )
    raise ValueError(f"Unknown Fig. 2D series {series!r}")


# ---------------------------------------------------------------------------
# Fig. 3B (p. 4): switching likelihood vs l_g (sites)
# ---------------------------------------------------------------------------

def read_fig3b(img: np.ndarray, at: Sequence[float]) -> list[dict]:
    return C.read_fig3b_circles(img, at)


# ---------------------------------------------------------------------------
# Fig. 3C Left (p. 4): domain length vs [fascin]/[alpha-actinin]
# ---------------------------------------------------------------------------

def read_fig3c_left_kmc(img: np.ndarray, at: Sequence[float], series: str) -> list[dict]:
    xa, ya, box, _ = C.calibrate_mpl_panel(
        img,
        x_tick_vals=[0.0, 1.0, 2.0],
        y_tick_vals=[10.0, 5.0, 0.0],
        x_scale="linear",
        y_scale="linear",
        y_side="left",
    )
    assert xa is not None
    c0, r0, c1, r1 = box
    r_g, g_g, b_g = C._rgb(img)
    if series == "fascin":
        m = (r_g > 200) & (g_g < 65) & (b_g < 65)
    elif series == "alpha_actinin":
        m = (r_g < 65) & (g_g > 200) & (b_g > 200)
    else:
        raise ValueError(f"Unknown Fig. 3C Left KMC series {series!r}")

    sub = m[r0:r1, c0:c1]
    dt = ndimage.distance_transform_edt(sub)
    rad_est = float(np.clip(dt.max(), 45.0, 52.0))

    xs = [float(x) for x in at] if at else [0.3, 0.5, 1.0, 1.5, 2.0]
    out: list[dict] = []
    for x in xs:
        col = xa.to_pixel(x)
        ci = int(round(col)) - c0
        sl = dt[:, max(0, ci - 6) : min(sub.shape[1], ci + 7)]
        mx = float(sl.max())
        if mx >= 0.7 * rad_est:
            rr = np.nonzero(sl >= mx - 0.5)[0] + r0
            rc = float(rr.mean())
            note = f"distance-transform circle center (R={mx:.1f} px)"
            half_px = 3.0
        else:
            rr_vis = np.nonzero(sub[:, max(0, ci - 15) : min(sub.shape[1], ci + 16)].any(axis=1))[0] + r0
            if len(rr_vis) == 0:
                out.append({"x": x, "y": None, "note": "occluded"})
                continue
            rc = float(rr_vis.min()) + rad_est
            note = f"partially occluded by cyan; top arc + R ({rad_est:.1f} px)"
            half_px = 6.0
        out.append(C._make_row(x=x, y=ya.to_value(rc), col=col, row=rc, ya=ya, half_px=half_px, note=note))
    return out


def _fast_square_opening(mask: np.ndarray, size: int) -> np.ndarray:
    """Separable O(HW) binary opening by a size x size square."""
    u = mask.astype(np.uint8)
    er = ndimage.minimum_filter1d(ndimage.minimum_filter1d(u, size=size, axis=0), size=size, axis=1)
    op = ndimage.maximum_filter1d(ndimage.maximum_filter1d(er, size=size, axis=0), size=size, axis=1)
    return op.astype(bool)


def read_fig3c_left_expt(img: np.ndarray, at: Sequence[float]) -> list[dict]:
    """Read cleanly visible black experimental squares + error bars in Fig. 3C Left."""
    xa, ya, box, _ = C.calibrate_mpl_panel(
        img,
        x_tick_vals=[0.0, 1.0, 2.0],
        y_tick_vals=[10.0, 5.0, 0.0],
        x_scale="linear",
        y_scale="linear",
        y_side="left",
    )
    assert xa is not None
    c0, r0, c1, r1 = box
    r_g, g_g, b_g = C._rgb(img)
    blk = (r_g < 45) & (g_g < 45) & (b_g < 45)
    valid = np.zeros_like(blk)
    valid[r0:r1, c0:c1] = blk[r0:r1, c0:c1]

    op = _fast_square_opening(valid[r0:r1, c0:c1], size=40)
    lbl, n = ndimage.label(op)
    raw_sq = []
    for k in range(1, n + 1):
        rr, cc = np.nonzero(lbl == k)
        rr = rr + r0
        cc = cc + c0
        w = cc.max() - cc.min() + 1
        h = rr.max() - rr.min() + 1
        if w > 120:
            c_left = float(cc.min() + 50.0)
            r_left = float(rr[cc < cc.min() + 95].mean())
            raw_sq.append((c_left, r_left, "leftmost bottom-left square"))
        elif w >= 55 and h >= 55:
            raw_sq.append((float(cc.mean()), float(rr.mean()), "unoccluded black square"))
    raw_sq.sort()

    col_col = ((r_g > 180) & (g_g < 80) & (b_g < 80)) | ((r_g < 80) & (g_g > 180) & (b_g > 180))
    out_pts: list[tuple[float, float, float, float, str]] = []
    for cm, rm, note in raw_sq:
        ci = int(round(cm))
        col_blk = valid[:, max(c0, ci - 5) : min(c1, ci + 6)].mean(axis=1) > 0.65
        col_bridge = col_col[:, max(c0, ci - 5) : min(c1, ci + 6)].any(axis=1)

        ru = max(r0, int(round(rm)) - 50)
        while ru > r0:
            if col_blk[ru - 1]:
                ru -= 1
            elif col_bridge[ru - 1]:
                look = ru - 1
                while look > r0 and col_bridge[look] and (ru - look) < 95:
                    look -= 1
                if look > r0 and col_blk[look]:
                    ru = look
                else:
                    break
            else:
                break
        rd = min(r1 - 1, int(round(rm)) + 50)
        while rd < r1 - 1:
            if col_blk[rd + 1]:
                rd += 1
            elif col_bridge[rd + 1]:
                look = rd + 1
                while look < r1 - 1 and col_bridge[look] and (look - rd) < 95:
                    look += 1
                if look < r1 - 1 and col_blk[look]:
                    rd = look
                else:
                    break
            else:
                break

        up_span = rm - ru
        dn_span = rd - rm
        if rd >= r1 - 4 and up_span > dn_span + 10:
            rd = rm + up_span
            note += "; lower error bar clipped at y=0, mirrored from upper"
        elif up_span > dn_span + 20 and rd > rm + 55:
            # Lower stem interrupted by crossing KMC line (e.g. x=1.05): mirror symmetrically
            rd = rm + up_span
            note += "; lower stem crossed by KMC line, mirrored symmetrically"
        out_pts.append((cm, rm, float(ru), float(rd), note))

    # Partially crossed square at x ~ 1.25 (cyan line crosses its center)
    c_s0, c_s1 = int(round(xa.to_pixel(1.18))), int(round(xa.to_pixel(1.32)))
    # Locate the clean vertical black error-bar stem below the cyan line (y in 2.0..3.1)
    r_stem0, r_stem1 = int(round(ya.to_pixel(3.1))), int(round(ya.to_pixel(2.0)))
    stem_prof = valid[r_stem0:r_stem1, c_s0:c_s1].sum(axis=0)
    if stem_prof.max() > 30:
        best_cs = np.nonzero(stem_prof >= stem_prof.max() * 0.9)[0] + c_s0
        cm_125 = float(best_cs.mean())
        ci = int(round(cm_125))
        col_blk = valid[:, ci - 5 : ci + 6].mean(axis=1) > 0.65
        col_br = col_col[:, ci - 5 : ci + 6].any(axis=1)
        # Find top and bottom of error bar along exact stem column ci, bridging across cyan line
        r_mid_s = int(round(ya.to_pixel(3.5)))
        # Walk up from above cyan
        above = np.nonzero(col_blk[int(round(ya.to_pixel(5.2))) : r_mid_s])[0]
        below = np.nonzero(col_blk[r_mid_s : int(round(ya.to_pixel(1.8)))])[0]
        if len(above) and len(below):
            ru = float(int(round(ya.to_pixel(5.2))) + above.min())
            rd = float(r_mid_s + below.max())
            rm_125 = 0.5 * (ru + rd)
            # Replace any partial top-half detection near 1.25
            out_pts = [p for p in out_pts if abs(p[0] - cm_125) > 45]
            out_pts.append((cm_125, rm_125, ru, rd, "black square crossed by cyan line; exact stem & symmetric error bar"))

    # Partially occluded square at x ~ 0.52 (red circle covers top-left)
    c_p0, c_p1 = int(round(xa.to_pixel(0.45))), int(round(xa.to_pixel(0.58)))
    r_off = int(round(ya.to_pixel(2.8)))
    sub_p = valid[r_off:r1, c_p0:c_p1]
    wide_p = np.nonzero(sub_p.sum(axis=1) >= 65)[0] + r_off
    if len(wide_p) >= 8:
        r_bot_p = float(wide_p.max())
        rm_p = r_bot_p - 50.0
        cc_p = np.nonzero(sub_p[wide_p - r_off].any(axis=0))[0] + c_p0
        cm_p = 0.5 * (float(cc_p.min()) + float(cc_p.max()))
        ci = int(round(cm_p))
        col_blk = valid[:, ci - 5 : ci + 6].mean(axis=1) > 0.65
        r_hi_s = int(round(ya.to_pixel(3.0)))
        ru_candidates = np.nonzero(col_blk[r_hi_s : int(round(rm_p)) - 45])[0]
        ru = float(r_hi_s + ru_candidates.min()) if len(ru_candidates) else rm_p - 50.0
        rd = rm_p + (rm_p - ru)
        if all(abs(cm_p - p[0]) > 40 for p in out_pts):
            out_pts.append((cm_p, rm_p, ru, rd, "black square partially behind red circle; symmetric error bar"))

    out_pts.sort(key=lambda t: t[0])
    return [
        C._make_row(
            x=xa.to_value(cm),
            y=ya.to_value(rm),
            col=cm,
            row=rm,
            ya=ya,
            half_px=4.0,
            bar_lo_row=rd,
            bar_hi_row=ru,
            bar_pad_px=3.0,
            note=note,
        )
        for cm, rm, ru, rd, note in out_pts
    ]


# ---------------------------------------------------------------------------
# Fig. 3C Right (p. 4): domain length vs mu_grow
# ---------------------------------------------------------------------------

def read_fig3c_right(img: np.ndarray, at: Sequence[float], series: str) -> list[dict]:
    xa, ya, box, _ = C.calibrate_mpl_panel(
        img,
        x_tick_vals=[0.0, 50.0, 100.0],
        y_tick_vals=[10.0, 5.0, 0.0],
        x_scale="linear",
        y_scale="linear",
        y_side="left",
    )
    assert xa is not None
    c0, r0, c1, r1 = box
    r_g, g_g, b_g = C._rgb(img)
    leg_r1 = int(round(ya.to_pixel(6.4)))
    leg_c0 = int(round(xa.to_pixel(40.0)))

    cyan = ((r_g < 65) & (g_g > 200) & (b_g > 200))[r0:r1, c0:c1].copy()
    red = ((r_g > 200) & (g_g < 65) & (b_g < 65))[r0:r1, c0:c1].copy()
    for m in (cyan, red):
        m[: max(0, leg_r1 - r0), max(0, leg_c0 - c0) :] = False

    dt_c = ndimage.distance_transform_edt(cyan)
    rad_est = float(np.clip(dt_c.max(), 47.0, 52.0))

    xs = [float(x) for x in at] if at else [20.0, 50.0, 100.0]
    out: list[dict] = []
    for x in xs:
        col = xa.to_pixel(x)
        ci = int(round(col)) - c0
        sl = dt_c[:, max(0, ci - 6) : min(cyan.shape[1], ci + 7)]
        mx = float(sl.max())
        rc_c = float((np.nonzero(sl >= mx - 0.5)[0] + r0).mean())
        if series == "alpha_actinin":
            rc = rc_c
            note = f"cyan circle center (R={mx:.1f} px)"
            half_px = 3.0
        elif series == "fascin":
            # Use narrow center strip (ci-3..ci+4) so steep connecting line from x=10 is not picked up
            rr_r = np.nonzero(red[:, max(0, ci - 3) : min(red.shape[1], ci + 4)].any(axis=1))[0] + r0
            if len(rr_r) and rr_r.min() < rc_c - rad_est + 5:
                rc = float(rr_r.min()) + rad_est
                note = f"red circle top arc peeking above cyan + R ({rad_est:.1f} px)"
                half_px = 5.0
            else:
                rc = rc_c
                note = "coincident with cyan circle center"
                half_px = 5.0
        else:
            raise ValueError(f"Unknown Fig. 3C Right series {series!r}")
        out.append(C._make_row(x=x, y=ya.to_value(rc), col=col, row=rc, ya=ya, half_px=half_px, note=note))
    return out


# ---------------------------------------------------------------------------
# Fig. 4B (p. 4) & Fig. S4C (p. 13): TIRF measured bar charts
# ---------------------------------------------------------------------------

def read_fig4b(img: np.ndarray, at: Sequence[float], series: str) -> list[dict]:
    color = "red" if series == "fascin" else "cyan"
    return C.read_bar_chart(img, y_tick_vals=[10.0, 5.0, 0.0], color=color, x_values=[1.5, 0.75])


def read_figS4c(img: np.ndarray, at: Sequence[float], series: str) -> list[dict]:
    color = "red" if series == "fascin" else "cyan"
    return C.read_bar_chart(img, y_tick_vals=[4.0, 3.0, 2.0, 1.0, 0.0], color=color, x_values=[3.0, 1.5])


# ---------------------------------------------------------------------------
# Fig. 5B (p. 5): short cross-linker domain length vs rho_short/rho_long
# ---------------------------------------------------------------------------

def read_fig5b(img: np.ndarray, at: Sequence[float], series: str) -> list[dict]:
    xa, ya, box, _ = C.calibrate_mathematica_panel(
        img,
        x_tick_vals=[0.1, 0.5, 1.0],
        y_tick_vals=[10.0, 8.0, 6.0, 4.0, 2.0, 0.0],
        x_scale="log10",
        y_scale="linear",
    )
    c0, r0, c1, r1 = box
    leg_r1 = int(round(ya.to_pixel(7.65)))
    leg_c1 = int(round(xa.to_pixel(1.35)))
    legend_rects = [(r0, leg_r1, c0, leg_c1)]

    if series == "dl27":
        mask = C.red_mask(img, r_min=150.0, diff=50.0)
        y_wins = {
            0.1: (-0.2, 0.65),
            0.2: (0.20, 0.405),
            0.5: (0.35, 0.68),
            1.0: (0.50, 1.35),
            1.5: (1.00, 2.20),
            2.0: (1.70, 3.20),
            3.0: (4.00, 5.65),
        }
        x_positions = [(v, v) for v in (0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0)]
        return C.read_mathematica_series_at_x(
            img, xa, ya, box, mask, x_positions,
            legend_rects=legend_rects, y_windows=y_wins, cap_half_w=34, filled_rad=14,
        )

    elif series == "dl100":
        mask = C.red_mask(img, r_min=150.0, diff=50.0)
        y_wins = {
            0.1: (-0.2, 0.65),
            0.2: (0.405, 0.72),
            0.5: (0.72, 1.15),
            1.0: (1.45, 2.35),
            1.5: (3.20, 5.50),
            2.0: (3.80, 7.50),
            3.0: (5.70, 10.90),
        }
        x_positions = [(v, v) for v in (0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0)]
        return C.read_mathematica_series_at_x(
            img, xa, ya, box, mask, x_positions, legend_rects=legend_rects, y_windows=y_wins, cap_half_w=34
        )

    elif series == "expt":
        H, W = img.shape[:2]
        blk = np.zeros((H, W), dtype=bool)
        raw_blk = C.black_mask(img, thresh=95.0, sat_max=35.0)
        blk[r0 + 20 : r1 - 6, c0:c1] = raw_blk[r0 + 20 : r1 - 6, c0:c1]
        for lr0, lr1, lc0, lc1 in legend_rects:
            blk[lr0:lr1, lc0:lc1] = False
        blk[r0 + 20 : leg_r1, int(round(xa.to_pixel(0.95))) : c1] = raw_blk[
            r0 + 20 : leg_r1, int(round(xa.to_pixel(0.95))) : c1
        ]
        blk[r0 : int(round(ya.to_pixel(8.95))), c0 : int(round(xa.to_pixel(1.25)))] = False

        # Horizontal opening on 2D strip (2x22) cleanly separates horizontal caps/square edges from diagonal lines
        h_blk = ndimage.binary_opening(blk, structure=np.ones((2, 22), bool))
        lbl, n = ndimage.label(h_blk)
        segs = []
        for k in range(1, n + 1):
            rr, cc = np.nonzero(lbl == k)
            w = int(cc.max() - cc.min() + 1)
            if 24 <= w <= 95:
                segs.append((float(cc.mean()), float(rr.mean()), w))
        segs.sort()
        groups: list[list[tuple[float, float, int]]] = []
        for s in segs:
            if not groups or abs(s[0] - np.mean([g[0] for g in groups[-1]])) > 22:
                groups.append([s])
            else:
                groups[-1].append(s)

        out: list[dict] = []
        for grp in groups:
            if len(grp) < 2:
                continue
            cm = float(np.mean([g[0] for g in grp]))
            all_r = sorted(g[1] for g in grp)
            caps = sorted(g[1] for g in grp if g[2] >= 60)
            sq_edges = sorted(g[1] for g in grp if 40 <= g[2] < 60)

            if all_r[-1] - all_r[0] <= 44 and abs(ya.to_value(0.5 * (all_r[0] + all_r[-1]))) < 0.15:
                center_r = 0.5 * (all_r[0] + all_r[-1])
                top_r = bot_r = center_r
                note = "open square at y~0; zero/hidden error bar"
            elif len(caps) >= 2:
                top_r, bot_r = caps[0], caps[-1]
                center_r = 0.5 * (top_r + bot_r)
                note = f"symmetric black error-bar caps (span={bot_r - top_r:.1f} px)"
            elif len(sq_edges) >= 2:
                center_r = 0.5 * (sq_edges[0] + sq_edges[-1])
                if len(caps) == 1:
                    half_eb = abs(caps[0] - center_r)
                    top_r, bot_r = center_r - half_eb, center_r + half_eb
                    note = "upper cap touches legend; mirrored symmetrically from lower cap & square center"
                else:
                    top_r, bot_r = all_r[0], all_r[-1]
                    note = "open square edges"
            else:
                top_r, bot_r = all_r[0], all_r[-1]
                center_r = 0.5 * (top_r + bot_r)
                note = "outermost black horizontal segments"

            out.append(
                C._make_row(
                    x=xa.to_value(cm),
                    y=ya.to_value(center_r),
                    col=cm,
                    row=center_r,
                    ya=ya,
                    half_px=3.0,
                    bar_lo_row=bot_r,
                    bar_hi_row=top_r,
                    bar_pad_px=2.5,
                    note=note,
                )
            )
        return out

    raise ValueError(f"Unknown Fig. 5B series {series!r}")


# ---------------------------------------------------------------------------
# Fig. 5C (p. 5): AFINES domain length vs k_grow (short & long)
# ---------------------------------------------------------------------------

def read_fig5c(img: np.ndarray, at: Sequence[float], series: str) -> list[dict]:
    xa, ya, box, _ = C.calibrate_mathematica_panel(
        img,
        x_tick_vals=[5.0, 10.0, 50.0, 100.0],
        y_tick_vals=[5.0, 4.0, 3.0, 2.0, 1.0],
        x_scale="log10",
        y_scale="linear",
        x_side="top",
    )
    c0, r0, c1, r1 = box
    leg_r1 = int(round(ya.to_pixel(4.25)))
    leg_c0 = int(round(xa.to_pixel(5.5)))
    leg_c1 = int(round(xa.to_pixel(20.0)))
    legend_rects = [(r0, leg_r1, leg_c0, leg_c1)]

    x_positions = [(v, v) for v in (2.5, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0)]
    if series == "short":
        mask = C.red_mask(img, r_min=150.0, diff=50.0)
    elif series == "long":
        mask = C.cyan_mask(img, gb_min=125.0, diff=35.0)
    else:
        raise ValueError(f"Unknown Fig. 5C series {series!r}")
    return C.read_mathematica_series_at_x(
        img, xa, ya, box, mask, x_positions, legend_rects=legend_rects, cap_half_w=34
    )


# ---------------------------------------------------------------------------
# Fig. 5E (p. 5): AFINES domain length vs k_grow at off1_ratio = 2 (top) & 4 (bot)
# ---------------------------------------------------------------------------

def read_fig5e(img: np.ndarray, at: Sequence[float], panel: str, series: str) -> list[dict]:
    panel_idx = 0 if panel == "top" else 1
    xa, ya, box, _ = C.calibrate_mathematica_panel(
        img,
        x_tick_vals=[10.0, 20.0, 50.0, 100.0],
        y_tick_vals=[4.0, 3.0, 2.0, 1.0],
        x_scale="log10",
        y_scale="linear",
        x_side="bottom",
        y_side="left",
        panel_idx=panel_idx,
    )
    c0, r0, c1, r1 = box
    legend_rects = []
    if panel == "top":
        legend_rects.append(
            (r0, int(round(ya.to_pixel(3.25))), int(round(xa.to_pixel(11.0))), int(round(xa.to_pixel(25.0))))
        )

    x_positions = [(v, v) for v in (10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0)]
    if series == "short":
        mask = C.red_mask(img, r_min=150.0, diff=50.0)
    elif series == "long":
        mask = C.cyan_mask(img, gb_min=125.0, diff=35.0)
    else:
        raise ValueError(f"Unknown Fig. 5E series {series!r}")
    return C.read_mathematica_series_at_x(
        img, xa, ya, box, mask, x_positions, legend_rects=legend_rects, cap_half_w=34
    )


# ---------------------------------------------------------------------------
# Fig. S2 (p. 11): AFINES calibration (domain length vs len diff, shifted +0.01 um)
# ---------------------------------------------------------------------------

def read_figS2(img: np.ndarray, at: Sequence[float], series: str) -> list[dict]:
    xa, ya, box, _ = C.calibrate_mathematica_panel(
        img,
        x_tick_vals=[0.01, 0.05, 0.10, 0.50, 1.0],
        y_tick_vals=[6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0],
        x_scale="log10",
        y_scale="linear",
    )
    x_positions = [(v, v) for v in (0.01, 0.037, 0.06, 0.11, 0.31, 0.51, 1.01)]
    if series == "short":
        mask = C.red_mask(img, r_min=150.0, diff=50.0)
    elif series == "long":
        mask = C.cyan_mask(img, gb_min=125.0, diff=35.0)
    else:
        raise ValueError(f"Unknown Fig. S2 series {series!r}")
    return C.read_mathematica_series_at_x(img, xa, ya, box, mask, x_positions, cap_half_w=34)


# ---------------------------------------------------------------------------
# Fig. S3A Left & Right, Fig. S3B (p. 12): KMC histogram curves
# ---------------------------------------------------------------------------

def read_figS3a_left(img: np.ndarray, at: Sequence[float], series: str) -> list[dict]:
    xa, ya, box, _ = C.calibrate_mpl_panel(
        img,
        x_tick_vals=[0.0, 5.0, 10.0, 15.0],
        y_tick_vals=[0.3, 0.2, 0.1, 0.0],
        x_scale="linear",
        y_scale="linear",
        y_side="left",
    )
    assert xa is not None
    c0, r0, c1, r1 = box
    legend_rects = [(r0, int(round(ya.to_pixel(0.22))), int(round(xa.to_pixel(3.8))), c1)]
    if series == "same":
        m = C.rgb_dist_mask(img, (31, 119, 180), tol=55.0)
        xs = [float(x) for x in at] if at else [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0]
    elif series == "diff":
        m = C.rgb_dist_mask(img, (255, 127, 14), tol=55.0)
        xs = [float(x) for x in at] if at else [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0]
    else:
        raise ValueError(f"Unknown Fig. S3A Left series {series!r}")
    return C.read_mpl_curve_at(img, xa, ya, box, m, xs, legend_rects=legend_rects)


def read_figS3a_right(img: np.ndarray, at: Sequence[float], series: str) -> list[dict]:
    xa, ya, box, _ = C.calibrate_mpl_panel(
        img,
        x_tick_vals=[0.0, 5.0, 10.0, 15.0, 20.0],
        y_tick_vals=[0.2, 0.1, 0.0],
        x_scale="linear",
        y_scale="linear",
        y_side="right",
    )
    assert xa is not None
    if series == "same":
        m = C.rgb_dist_mask(img, (31, 119, 180), tol=55.0)
        xs = [float(x) for x in at] if at else [4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    elif series == "diff":
        m = C.rgb_dist_mask(img, (255, 127, 14), tol=55.0)
        xs = [float(x) for x in at] if at else [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 14.0, 16.0, 18.0]
    else:
        raise ValueError(f"Unknown Fig. S3A Right series {series!r}")
    return C.read_mpl_curve_at(img, xa, ya, box, m, xs)


def read_figS3b(img: np.ndarray, at: Sequence[float]) -> list[dict]:
    xa, ya, box, _ = C.calibrate_mpl_panel(
        img,
        x_tick_vals=[0.0, 50.0, 100.0, 150.0, 200.0, 250.0],
        y_tick_vals=[0.02, 0.01, 0.0],
        x_scale="linear",
        y_scale="linear",
        y_side="left",
    )
    assert xa is not None
    c0, r0, c1, r1 = box
    legend_rects = [(r0, int(round(ya.to_pixel(0.008))), int(round(xa.to_pixel(155.0))), c1)]
    m = C.cyan_mask(img, gb_min=190.0, diff=100.0)
    xs = [float(x) for x in at] if at else [
        5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 245.0
    ]
    return C.read_mpl_curve_at(img, xa, ya, box, m, xs, legend_rects=legend_rects)


# ---------------------------------------------------------------------------
# Fig. S3C (p. 12): KMC master collapse eps_predicted vs mu_grow(sites/s)/mu_on
# ---------------------------------------------------------------------------

_FIGS3C_SERIES = {
    "5nms":   ((31, 119, 180),  [0.7143, 1.4286]),
    "10nms":  ((255, 127, 14),  [1.4286]),
    "20nms":  ((44, 160, 44),   [1.1429, 2.8571]),
    "30nms":  ((214, 39, 40),   [0.8571, 1.7143, 4.2857, 8.5714, 17.1429, 42.8571]),
    "40nms":  ((148, 103, 189), [2.2857, 5.7143, 11.4286, 22.8571]),
    "50nms":  ((140, 86, 75),   [1.4286, 7.1429, 71.4286]),
    "100nms": ((227, 119, 194), [2.8571, 5.7143, 14.2857, 28.5714, 57.1429, 142.8571]),
}


def read_figS3c(img: np.ndarray, at: Sequence[float], series: str) -> list[dict]:
    xa, ya, box, _ = C.calibrate_mpl_panel(
        img,
        x_tick_vals=[1.0, 10.0, 100.0, 1000.0],
        y_tick_vals=[8.0, 6.0, 4.0, 2.0, 0.0],
        x_scale="log10",
        y_scale="linear",
        y_side="left",
    )
    assert xa is not None
    c0, r0, c1, r1 = box
    rgb, default_xs = _FIGS3C_SERIES[series]
    sub_img = img[r0:r1, c0:c1]
    sub = C.rgb_dist_mask(sub_img, rgb, tol=36.0)
    # Mask legend box (upper right: x > 14, y > 3.6)
    leg_r1 = max(0, int(round(ya.to_pixel(3.6))) - r0)
    leg_c0 = max(0, int(round(xa.to_pixel(14.0))) - c0)
    sub[:leg_r1, leg_c0:] = False

    dt = ndimage.distance_transform_edt(sub)
    R_full = 38.0

    xs_hint = [float(x) for x in at] if at else default_xs
    out: list[dict] = []
    for x_hint in xs_hint:
        c_est = int(round(xa.to_pixel(x_hint)))
        ci = c_est - c0
        win_c0, win_c1 = max(0, ci - 24), min(sub.shape[1], ci + 25)
        dt_sub = dt[:, win_c0:win_c1]
        mx = float(dt_sub.max())
        if mx >= 26.0:
            rr_pk, cc_pk = np.nonzero(dt_sub >= mx - 0.8)
            rc = r0 + float(rr_pk.mean())
            cc = c0 + win_c0 + float(cc_pk.mean())
            x_meas = xa.to_value(cc)
            note = f"circle center (dt_max={mx:.1f} px)"
            half_px = 3.0
        else:
            strip = sub[:, max(0, ci - 8) : min(sub.shape[1], ci + 9)]
            rr_vis = np.nonzero(strip.any(axis=1))[0]
            if len(rr_vis) == 0:
                out.append({"x": x_hint, "y": None, "note": "fully occluded by another series"})
                continue
            runs = C._runs_1d(strip.any(axis=1), max_gap=4)
            best = max(runs, key=lambda t: t[1] - t[0])
            top_e, bot_e = r0 + best[0], r0 + best[1]
            span = bot_e - top_e
            if span >= 55:
                rc = 0.5 * (top_e + bot_e)
                half_px = 4.0
                note = f"partially occluded circle, full vertical span ({span} px)"
            else:
                above_white = img[max(r0, top_e - 6), c_est].mean() > 235
                below_white = img[min(r1 - 1, bot_e + 6), c_est].mean() > 235
                if above_white and not below_white:
                    rc = float(top_e) + R_full
                    note = f"upper crescent + R ({R_full:.0f} px)"
                    half_px = 6.0
                elif below_white and not above_white:
                    rc = float(bot_e) - R_full
                    note = f"lower crescent - R ({R_full:.0f} px)"
                    half_px = 6.0
                else:
                    rc = 0.5 * (top_e + bot_e)
                    note = f"partial marker midpoint (span={span} px)"
                    half_px = 6.0
            cc = float(c_est)
            x_meas = xa.to_value(cc)
        out.append(C._make_row(x=x_meas, y=ya.to_value(rc), col=cc, row=rc, ya=ya, half_px=half_px, note=note))
    return out


# ---------------------------------------------------------------------------
# Fig. S5 (p. 14): AFINES domain length vs k_grow (1500-1600 s window)
# ---------------------------------------------------------------------------

def read_figS5(img: np.ndarray, at: Sequence[float], series: str) -> list[dict]:
    xa, ya, box, _ = C.calibrate_mathematica_panel(
        img,
        x_tick_vals=[0.01, 0.02, 0.05, 0.10],
        y_tick_vals=[3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.0],
        x_scale="log10",
        y_scale="linear",
    )
    x_positions = [(v, v) for v in (0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.10)]
    if series == "short":
        mask = C.red_mask(img, r_min=150.0, diff=50.0)
    elif series == "long":
        mask = C.cyan_mask(img, gb_min=125.0, diff=35.0)
    else:
        raise ValueError(f"Unknown Fig. S5 series {series!r}")
    return C.read_mathematica_series_at_x(img, xa, ya, box, mask, x_positions, cap_half_w=12)
