"""Calibrated readers for Fig. 4a-c, Fig. 6b-d, Supp Fig. 6, and Supp Fig. 9 (medians, boxplots, scatter)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy.ndimage import maximum_filter
from scipy.signal import fftconvolve

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as C


def read_fig4a(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    """Read Fig. 4a horizontal median bars for Domain diameter (um) at c_CP = 60, 120 nM."""
    dark_all, gray, cl, rt, cr, rb = C.find_box(image, thresh=150.0)
    sat = image.max(axis=2) - image.min(axis=2)
    ndark = dark_all & (sat < 20)
    yt = C.find_inward_ticks(ndark, cl, rt, cr, rb, "y_left", band=(4, 18), margin=12)
    # 4 interior ticks: 8, 6, 4, 2 um (tick 10 sits just below top spine and merges slightly)
    yt_m = C.match_ticks_regular(yt, 4, tol=6.0)
    ya = C.verify_axis("linear", list(zip(yt_m, [8.0, 6.0, 4.0, 2.0])), max_residual_px=2.5)
    xt = C.find_inward_ticks(ndark, cl, rt, cr, rb, "x_top", band=(4, 18), margin=12)
    if len(xt) != 2:
        raise RuntimeError(f"Expected 2 x-ticks in Fig. 4a, found {xt}")

    if series == "myosin":
        mask = (image[:, :, 0] > 200) & (image[:, :, 1] > 180) & (image[:, :, 2] < 110)
    else:
        mask = (sat < 15) & (gray > 180) & (gray < 380)

    out = []
    for cp, xc in zip([60.0, 120.0], xt):
        wings = list(range(int(xc - 420), int(xc - 280))) + list(range(int(xc + 280), int(xc + 420)))
        prof = mask[int(rt[2]) + 20 : int(rb[1]) - 20, wings].mean(axis=1)
        bars = C.clusters(np.where(prof > 0.6)[0] + int(rt[2]) + 20, gap=3)
        if not bars:
            raise RuntimeError(f"No {series} median bar found in Fig. 4a at {cp} nM")
        r_mid, r0, r1 = bars[0]
        thick = float(r1 - r0 + 1)
        out.append(
            C.make_point_reading(
                cp,
                xc,
                r_mid,
                ya,
                half_px=max(0.5 * thick, 2.0) + ya.residual_px,
                note=f"Fig4a {series} median bar, thick={thick:.0f}px, cal_res={ya.residual_px:.2f}px",
            )
        )
    return out


def read_fig4b(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    """Read Fig. 4b horizontal median bars for Curvature (1/um, negative) at c_CP = 60, 120 nM."""
    dark_all, gray, cl, rt, cr, rb = C.find_box(image, thresh=150.0)
    sat = image.max(axis=2) - image.min(axis=2)
    ndark = dark_all & (sat < 20)
    yt = C.find_inward_ticks(ndark, cl, rt, cr, rb, "y_left", band=(4, 18), margin=12)
    yt_m = C.match_ticks_regular(yt + [rb[0]], 4, tol=6.0)
    ya = C.verify_axis("linear", list(zip(yt_m, [-0.6, -0.4, -0.2, 0.0])), max_residual_px=2.5)
    xt = C.find_inward_ticks(ndark, cl, rt, cr, rb, "x_top", band=(4, 18), margin=12)
    if len(xt) != 2:
        raise RuntimeError(f"Expected 2 x-ticks in Fig. 4b, found {xt}")

    if series == "myosin":
        mask = (image[:, :, 0] > 200) & (image[:, :, 1] > 180) & (image[:, :, 2] < 110)
    else:
        mask = (sat < 15) & (gray > 180) & (gray < 380)

    out = []
    for cp, xc in zip([60.0, 120.0], xt):
        wings = list(range(int(xc - 420), int(xc - 280))) + list(range(int(xc + 280), int(xc + 420)))
        prof = mask[int(rt[2]) + 20 : int(rb[1]) - 20, wings].mean(axis=1)
        bars = C.clusters(np.where(prof > 0.6)[0] + int(rt[2]) + 20, gap=3)
        if not bars:
            raise RuntimeError(f"No {series} median bar found in Fig. 4b at {cp} nM")
        r_mid, r0, r1 = bars[0]
        thick = float(r1 - r0 + 1)
        out.append(
            C.make_point_reading(
                cp,
                xc,
                r_mid,
                ya,
                half_px=max(0.5 * thick, 2.0) + ya.residual_px,
                note=f"Fig4b {series} median bar, thick={thick:.0f}px, cal_res={ya.residual_px:.2f}px",
            )
        )
    return out


def read_fig4c(image: np.ndarray, at: list[float], cp_nM: int) -> list[dict]:
    """Read Fig. 4c scatter markers (Curvature vs Domain diameter) via legend-template correlation."""
    dark_all, gray, cl, rt, cr, rb = C.find_box(image, thresh=150.0)
    sat = image.max(axis=2) - image.min(axis=2)
    ndark = dark_all & (sat < 20)
    yt = C.find_inward_ticks(ndark, cl, rt, cr, rb, "y_left", band=(4, 18), margin=12)
    xt = C.find_inward_ticks(ndark, cl, rt, cr, rb, "x_bot", band=(4, 18), margin=12)
    ya = C.verify_axis(
        "linear",
        list(zip(C.match_ticks_regular(yt + [rt[0], rb[0]], 6, tol=6.0), [-0.5, -0.4, -0.3, -0.2, -0.1, 0.0])),
        max_residual_px=2.5,
    )
    xa = C.verify_axis(
        "linear",
        list(zip(C.match_ticks_regular(xt + [cl[0]], 5, tol=6.0), [1.0, 3.0, 5.0, 7.0, 9.0])),
        max_residual_px=3.0,
    )

    g = (gray // 3).astype(np.int16)
    # Tight column band containing ONLY the legend marker icons (x in [5.95, 6.62], well left of '60 nM CP')
    lc0, lc1 = int(xa.to_pixel(5.95)), int(xa.to_pixel(6.62))
    if cp_nM == 60:
        mask = (sat < 15) & (g >= 172) & (g <= 228)
        lr0, lr1 = int(ya.to_pixel(-0.475)), int(ya.to_pixel(-0.415))
        thresh = 0.68
    else:
        mask = (sat < 15) & (g >= 50) & (g <= 125)
        lr0, lr1 = int(ya.to_pixel(-0.405)), int(ya.to_pixel(-0.335))
        thresh = 0.64

    sub = mask[lr0:lr1, lc0:lc1]
    rr, cc = np.where(sub)
    if len(rr) == 0:
        raise RuntimeError(f"Legend marker not found for Fig. 4c {cp_nM} nM")
    tmpl = sub[int(rr.min()) : int(rr.max()) + 1, int(cc.min()) : int(cc.max()) + 1].astype(np.float32)

    work = mask.copy().astype(np.float32)
    work[: int(rt[2]) + 6, :] = 0
    work[int(rb[1]) - 4 :, :] = 0
    work[:, : int(cl[2]) + 6] = 0
    work[:, int(cr[1]) - 6 :] = 0
    # Exclude legend area
    work[: int(ya.to_pixel(-0.32)), int(xa.to_pixel(5.7)) :] = 0

    work2 = work[::2, ::2]
    tmpl2 = tmpl[::2, ::2]
    t2_sum = float(tmpl2.sum())
    score = fftconvolve(work2, tmpl2[::-1, ::-1], mode="same") / max(t2_sum, 1.0)
    mx = maximum_filter(score, size=12)
    pks = np.argwhere((score == mx) & (score >= thresh))

    # Uncertainty: half-pixel + calibration residual + 4 px marker centroid tolerance
    half_px = 4.0 + ya.residual_px
    pts = []
    for r2, c2 in pks:
        r_f, c_f = float(r2 * 2 + 0.5), float(c2 * 2 + 0.5)
        x_val = float(xa.to_value(c_f))
        pts.append((x_val, c_f, r_f, float(score[r2, c2])))
    pts.sort(key=lambda t: (t[0], t[2]))

    out = []
    for x_val, c_f, r_f, sc in pts:
        out.append(
            C.make_point_reading(
                x_val,
                c_f,
                r_f,
                ya,
                half_px=half_px,
                note=f"Fig4c {cp_nM}nM marker (corr={sc:.2f}), cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
            )
        )
    return out


def _read_fig6_boxplot(
    image: np.ndarray,
    y_vals: list[float],
    panel_name: str,
    inset_masks: list[tuple[float, float, float, float]] | None = None,
) -> list[dict]:
    """Read Fig. 6b/c/d boxplots (red median + grey non-outlier whisker range) across 3 conditions."""
    dark_all, gray, cl, rt, cr, rb = C.find_box(image, thresh=150.0)
    sat = image.max(axis=2) - image.min(axis=2)
    ndark = dark_all & (sat < 20)
    yt = C.find_inward_ticks(ndark, cl, rt, cr, rb, "y_left", band=(4, 18), margin=12)
    xt = C.find_inward_ticks(ndark, cl, rt, cr, rb, "x_bot", band=(4, 18), margin=12)
    ya = C.verify_axis("linear", list(zip(C.match_ticks_regular(yt, len(y_vals), tol=6.0), y_vals)), max_residual_px=2.5)
    if len(xt) != 3:
        raise RuntimeError(f"Expected 3 x-ticks in {panel_name}, found {xt}")

    red = (image[:, :, 0] > 140) & (image[:, :, 1] < 70) & (image[:, :, 2] < 70)
    gry = (sat < 18) & (gray < 520)

    # Mask out micrograph insets (given in fractional box coords [fx0, fx1, y_lo, y_hi])
    if inset_masks:
        w_box = cr[0] - cl[0]
        for fx0, fx1, y_lo, y_hi in inset_masks:
            ic0, ic1 = int(cl[0] + fx0 * w_box), int(cl[0] + fx1 * w_box)
            ir0, ir1 = sorted((int(ya.to_pixel(y_hi)), int(ya.to_pixel(y_lo))))
            red[max(0, ir0) : ir1, max(0, ic0) : ic1] = False
            gry[max(0, ir0) : ir1, max(0, ic0) : ic1] = False

    cond_names = ["CP Dox", "no CP Dox", "CP no Dox"]
    r_top, r_bot = int(rt[2]) + 8, int(rb[1]) - 5
    out = []
    for idx, (xc, cname) in enumerate(zip(xt, cond_names), start=1):
        # 1. Red horizontal median line
        rprof = red[r_top:r_bot, int(xc - 135) : int(xc + 135)].mean(axis=1)
        rc = C.clusters(np.where(rprof > 0.30)[0] + r_top, gap=3)
        # 2. Grey box top/bottom edges (outer box wings xc +/- 90..135)
        box_cols = list(range(int(xc - 135), int(xc - 88))) + list(range(int(xc + 88), int(xc + 135)))
        bprof = gry[r_top:r_bot, box_cols].mean(axis=1)
        bc = C.clusters(np.where(bprof > 0.55)[0] + r_top, gap=3)
        # 3. Whisker caps + box edges in outer whisker-cap wings xc +/- 52..74 (avoids scatter circles near xc)
        w_cols = list(range(int(xc - 74), int(xc - 52))) + list(range(int(xc + 52), int(xc + 74)))
        wprof = gry[r_top:r_bot, w_cols].mean(axis=1)
        wc = [c for c in C.clusters(np.where(wprof > 0.55)[0] + r_top, gap=3) if (c[2] - c[1] + 1) <= 18]

        if rc:
            med_row = rc[0][0]
            med_thick = float(rc[0][2] - rc[0][1] + 1)
        elif bc:
            # Degenerate zero-height box (e.g. Fig. 6b CP Dox where median = Q1 = Q3 = 0 and red line is covered by dots)
            med_row = bc[0][0]
            med_thick = float(bc[0][2] - bc[0][1] + 1)
        else:
            raise RuntimeError(f"Could not find median for {panel_name} cond {idx}")

        all_rows = [c[0] for c in bc] + [c[0] for c in wc] + [med_row]
        row_hi = min(all_rows)  # highest value (smallest row index)
        row_lo = max(all_rows)  # lowest value (largest row index)
        q_vals = sorted(float(ya.to_value(c[0])) for c in bc)
        q_str = f"Q1={q_vals[0]:.3g}, Q3={q_vals[-1]:.3g}" if len(q_vals) >= 2 else "degenerate box"

        out.append(
            C.make_point_reading(
                float(idx),
                xc,
                med_row,
                ya,
                half_px=max(0.5 * med_thick, 2.0) + ya.residual_px,
                bar_row_lo=row_lo,
                bar_row_hi=row_hi,
                note=f"{panel_name} {cname} (x={idx}): {q_str}, whiskers=[{ya.to_value(row_lo):.3g}, {ya.to_value(row_hi):.3g}], cal_res={ya.residual_px:.2f}px",
            )
        )
    return out


def read_fig6b(image: np.ndarray, at: list[float]) -> list[dict]:
    return _read_fig6_boxplot(
        image,
        [20.0, 15.0, 10.0, 5.0, 0.0],
        "Fig6b",
        inset_masks=[(0.55, 0.80, 16.5, 24.5)],
    )


def read_fig6c(image: np.ndarray, at: list[float]) -> list[dict]:
    return _read_fig6_boxplot(
        image,
        [2.5, 2.0, 1.5, 1.0, 0.5],
        "Fig6c",
        inset_masks=[(0.05, 0.36, 1.75, 2.65), (0.50, 0.78, 0.45, 1.25)],
    )


def read_fig6d(image: np.ndarray, at: list[float]) -> list[dict]:
    return _read_fig6_boxplot(
        image,
        [0.5, 0.4, 0.3, 0.2, 0.1, 0.0],
        "Fig6d",
        inset_masks=[(0.05, 0.36, 0.01, 0.19)],
    )


def read_suppfig6(image: np.ndarray, at: list[float]) -> list[dict]:
    """Read Supp Fig. 6 domain-count boxplots for 1=protrusion, 2=concave domain, 3=flat domain."""
    dark_all, gray, cl, rt, cr, rb = C.find_box(image, thresh=170.0)
    sat = image.max(axis=2) - image.min(axis=2)
    ndark = dark_all & (sat < 20)
    yt = C.find_inward_ticks(ndark, cl, rt, cr, rb, "y_left", band=(4, 18), margin=12)
    xt = C.find_inward_ticks(ndark, cl, rt, cr, rb, "x_bot", band=(4, 18), margin=12)
    ya = C.verify_axis(
        "linear", list(zip(C.match_ticks_regular(yt, 5, tol=6.0), [80.0, 60.0, 40.0, 20.0, 0.0])), max_residual_px=2.5
    )
    if len(xt) != 3:
        raise RuntimeError(f"Expected 3 x-ticks in Supp Fig. 6, found {xt}")

    gry = (sat < 15) & (gray > 280) & (gray < 600)
    blk = (sat < 15) & (gray < 280)
    r_top, r_bot = int(rt[2]) + 15, int(rb[1]) - 15
    cats = ["protrusion", "concave domain", "flat domain"]

    out = []
    for idx, (xc, cname) in enumerate(zip(xt, cats), start=1):
        box_cols = list(range(int(xc - 170), int(xc - 80))) + list(range(int(xc + 80), int(xc + 170)))
        bprof = gry[r_top:r_bot, box_cols].mean(axis=1)
        bc = C.clusters(np.where(bprof > 0.6)[0] + r_top, gap=3)
        if len(bc) != 3:
            raise RuntimeError(f"Expected 3 horizontal box lines (Q3, median, Q1) in Supp Fig 6 cat {idx}, got {bc}")
        q3_row, med_row, q1_row = bc[0][0], bc[1][0], bc[2][0]
        med_thick = float(bc[1][2] - bc[1][1] + 1)

        w_cols = list(range(int(xc - 70), int(xc - 25))) + list(range(int(xc + 25), int(xc + 70)))
        wprof = blk[r_top:r_bot, w_cols].mean(axis=1)
        wc = C.clusters(np.where(wprof > 0.6)[0] + r_top, gap=3)
        all_rows = [c[0] for c in wc] + [q3_row, q1_row]
        row_hi, row_lo = min(all_rows), max(all_rows)
        out.append(
            C.make_point_reading(
                float(idx),
                xc,
                med_row,
                ya,
                half_px=max(0.5 * med_thick, 2.0) + ya.residual_px,
                bar_row_lo=row_lo,
                bar_row_hi=row_hi,
                note=f"SuppFig6 {cname} (x={idx}): Q1={ya.to_value(q1_row):.1f}, Q3={ya.to_value(q3_row):.1f}, whiskers=[{ya.to_value(row_lo):.1f}, {ya.to_value(row_hi):.1f}], cal_res={ya.residual_px:.2f}px",
            )
        )
    return out


def read_suppfig9(image: np.ndarray, at: list[float], kind: str) -> list[dict]:
    """Read Supp Fig. 9 horizontal median bar at 40 nM CP ('diameter' or 'curvature')."""
    y_vals = [12.0, 10.0, 8.0, 6.0, 4.0, 2.0, 0.0] if kind == "diameter" else [0.25, 0.20, 0.15, 0.10, 0.05, 0.0]
    dark_all, gray, cl, rt, cr, rb = C.find_box(image, thresh=160.0)
    sat = image.max(axis=2) - image.min(axis=2)
    ndark = dark_all & (sat < 20)
    yt = C.find_inward_ticks(ndark, cl, rt, cr, rb, "y_left", band=(4, 18), margin=12)
    ya = C.verify_axis(
        "linear", list(zip(C.match_ticks_regular(yt + [rt[0], rb[0]], len(y_vals), tol=8.0), y_vals)), max_residual_px=3.0
    )
    wings = list(range(int(cl[2]) + 30, int(cl[2]) + 280)) + list(range(int(cr[1]) - 280, int(cr[1]) - 30))
    gry = (sat < 15) & (gray > 280) & (gray < 640)
    prof = gry[int(rt[2]) + 15 : int(rb[1]) - 15, wings].mean(axis=1)
    ms = C.clusters(np.where(prof > 0.7)[0] + int(rt[2]) + 15, gap=3)
    if not ms:
        raise RuntimeError(f"Could not find median line in Supp Fig. 9 ({kind})")
    r_mid, r0, r1 = ms[0]
    thick = float(r1 - r0 + 1)
    xc = 0.5 * (cl[0] + cr[0])
    # Also find min/max of scatter circles inside center column region for note
    circ = gry[int(rt[2]) + 10 : int(rb[1]) - 10, int(cl[0] + 0.2 * (cr[0] - cl[0])) : int(cr[0] - 0.2 * (cr[0] - cl[0]))]
    c_rows = np.where(circ.any(axis=1))[0]
    span_note = ""
    if len(c_rows):
        y_max = ya.to_value(float(c_rows[0] + int(rt[2]) + 10))
        y_min = ya.to_value(float(c_rows[-1] + int(rt[2]) + 10))
        span_note = f", scatter_span=[{y_min:.3g}, {y_max:.3g}]"
    return [
        C.make_point_reading(
            40.0,
            xc,
            r_mid,
            ya,
            half_px=max(0.5 * thick, 2.0) + ya.residual_px,
            note=f"SuppFig9 40nM {kind} median bar, thick={thick:.0f}px{span_note}, cal_res={ya.residual_px:.2f}px",
        )
    ]
