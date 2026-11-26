"""Calibrated reader for Fig. S3C (simulated KMC domain-length distributions P(l), semi-log)."""

from __future__ import annotations

import numpy as np

import _calib as C


def _calibrate(image: np.ndarray):
    # Closed rectangular frame: top ~313.5, bottom ~2188.5, left ~560.0, right ~2922.0
    dark, gray, h_lines, v_lines = C.find_l_spines(
        image, thresh=140, row_frac=0.60, col_frac=0.60, row_range=(250, 2250), col_range=(450, 3000)
    )
    top_spine = min(h_lines, key=lambda t: t[0])
    bot_spine = max(h_lines, key=lambda t: t[0])
    left_spine = min(v_lines, key=lambda t: t[0])
    right_spine = max(v_lines, key=lambda t: t[0])

    # Bottom x-axis (l in um): left_spine=0, inward ticks at 2,4,6,8,[10],12,14, right_spine=16
    xt_in = C.find_ticks(dark, bot_spine, (int(left_spine[2]) + 20, int(right_spine[1]) - 20), "x", side=-1, band=(5, 22), frac=0.4, max_width=20)
    # Match against ideal 2..14 um grid interpolated between left_spine (0) and right_spine (16)
    x_marks = [(left_spine[0], 0.0), (right_spine[0], 16.0)]
    span_x = right_spine[0] - left_spine[0]
    for lv in (2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0):
        pred = left_spine[0] + (lv / 16.0) * span_x
        near = [t for t in xt_in if abs(t - pred) <= 15.0]
        if near:
            x_marks.append((min(near, key=lambda t: abs(t - pred)), lv))
    x_marks.sort(key=lambda t: t[1])
    x_axis = C.verify_axis("linear", x_marks, max_residual_px=4.5)

    # Y-axis (log10 P(l)): calibrate from printed decade labels 10^-1 .. 10^-5 AND inward decade ticks on right spine!
    # Right-inward ticks at 10^-2 (~934), 10^-3 (~1350), 10^-4 (~1773) + bottom spine 10^-5 (~2188.5)
    yt_r = C.find_ticks(dark, right_spine, (int(top_spine[2]) + 100, int(bot_spine[1]) - 50), "y", side=-1, band=(5, 22), frac=0.4, max_width=20)
    # Also printed decade label centers left of left_spine (10^-1, 10^-2, 10^-3, 10^-4, 10^-5)
    c_hi = int(left_spine[1]) - 20
    c_lo = max(0, c_hi - 340)
    trows = (gray[350 : int(bot_spine[2]) + 80, c_lo:c_hi] < 140).any(axis=1)
    lbl_cs = [c for c in C.clusters(np.where(trows)[0] + 350, gap=30) if (c[2] - c[1] + 1) >= 80]
    if len(lbl_cs) != 5:
        raise RuntimeError(f"Fig. S3C: expected 5 decade y-labels, found {lbl_cs}")
    # Anchor exact inward tick rows where detected (10^-3 and 10^-4 are completely unobstructed on right spine)
    # plus 10^-5 at bottom spine and check linearity against label centers
    y_marks = [(bot_spine[0], 1e-5)]
    step_lbl = (lbl_cs[-1][0] - lbl_cs[0][0]) / 4.0
    for k, val in enumerate([1e-1, 1e-2, 1e-3, 1e-4]):
        pred_r = bot_spine[0] - (4 - k) * step_lbl
        near = [t for t in yt_r if abs(t - pred_r) <= 18.0]
        if near:
            y_marks.append((min(near, key=lambda t: abs(t - pred_r)), val))
        else:
            # Fall back to exact decade offset from bottom spine using label vertical center shift
            y_marks.append((lbl_cs[k][0] + (bot_spine[0] - lbl_cs[-1][0]), val))
    y_marks.sort(key=lambda t: t[0])
    y_axis = C.verify_axis("log10", y_marks, max_residual_px=5.0)
    return x_axis, y_axis, (left_spine[0], top_spine[0], right_spine[0], bot_spine[0])


def read_series(image: np.ndarray, series: str, at: list[float]) -> list[dict]:
    """Read one of the 5 cleanly unoccluded Fig. S3C branches:
    fascin_r2p0 (top cyan), fascin_r1p0 (mid cyan), fascin_r0p5 (bot cyan),
    actinin_r1p0 (mid red), actinin_r2p0 (bot red).
    (Note: actinin_r0p5 top red branch is heavily occluded behind its symmetric twin fascin_r2p0.)
    """
    x_axis, y_axis, (c0, r0, c1, r1) = _calibrate(image)
    R = image[:, :, 0].astype(np.float32)
    G = image[:, :, 1].astype(np.float32)
    B = image[:, :, 2].astype(np.float32)

    is_cyan = series.startswith("fascin_")
    if is_cyan:
        mask = (B > 160) & (G > 140) & (R < 130)
    else:
        mask = (R > 160) & (R - G > 70) & (R - B > 70)

    # Mask outside frame and both legend boxes
    r_min, r_max = int(r0) + 12, int(r1) - 8
    c_min, c_max = int(c0) + 8, int(c1) - 8
    mask[:r_min, :] = False
    mask[r_max + 1 :, :] = False
    mask[:, :c_min] = False
    mask[:, c_max + 1 :] = False
    # Upper-right legend box (x > 3.2 um, P > 0.011) & lower-left legend box (x < 6.4 um, P < 1e-4)
    mask[int(r0) + 20 : int(round(y_axis.to_pixel(0.0105))), int(round(x_axis.to_pixel(3.0))) : c_max] = False
    mask[int(round(y_axis.to_pixel(1.1e-4))) : r_max, c_min : int(round(x_axis.to_pixel(6.5)))] = False

    # Expected approximate log10(P(l)) guide lines ONLY to disambiguate which of the 3 branches
    # (shallow / medium / steep) is picked at each x when multiple branches share the same color!
    # The actual read y-value comes 100% from the detected image pixels of that branch.
    # For the finite-lattice boundary spike at l > 14.3 um, pick by top/mid/bot vertical rank.
    guide_params = {
        "fascin_r2p0": (-2.42, -0.050),   # log10 P ~ a + b*l  (top cyan)
        "fascin_r1p0": (-2.05, -0.125),   # mid cyan
        "fascin_r0p5": (-2.00, -0.215),   # bot cyan
        "actinin_r1p0": (-2.05, -0.125),  # mid red
        "actinin_r2p0": (-2.00, -0.215),  # bot red
    }
    if series not in guide_params:
        raise ValueError(f"Unknown Fig. S3C series {series!r}")
    ga, gb = guide_params[series]

    out = []
    for xv in at:
        xv = float(xv)
        c_f = float(x_axis.to_pixel(xv))
        c = int(round(c_f))
        # Collect contiguous colored bands in a narrow window [c-8, c+8] (captures marker+line centers)
        lo_c, hi_c = max(c_min, c - 8), min(c_max, c + 9)
        prof = mask[r_min : r_max + 1, lo_c:hi_c].any(axis=1)
        segs = [(a + r_min, b + r_min) for a, b in C.runs_1d(prof, gap=10)]
        if not segs:
            out.append({"x": xv, "y": None, "found": False, "note": f"no {series} pixels at l={xv} um"})
            continue

        if xv >= 14.3:
            # Finite-lattice N=400 boundary spike column: order by vertical position (top=longest domain)
            segs.sort(key=lambda s: s[0])
            if series == "fascin_r2p0":
                seg = segs[0]
                # For the topmost cyan spike, take the actual marker peak center (top ~35 px of spike)
                sub_rr = np.where(mask[seg[0] : min(seg[1] + 1, seg[0] + 45), lo_c:hi_c])[0] + seg[0]
                r_est = float(sub_rr.mean()) + 0.5
                thick = 36.0
            elif series == "fascin_r1p0":
                # Middle spike sits around P ~ 1.6e-3 (row ~1260)
                pred_r = float(y_axis.to_pixel(1.6e-3))
                seg = min(segs, key=lambda s: abs(0.5 * (s[0] + s[1]) - pred_r))
                r_est = 0.5 * (seg[0] + seg[1]) + 0.5
                thick = float(seg[1] - seg[0] + 1)
            elif series == "fascin_r0p5":
                seg = segs[-1]
                r_est = 0.5 * (seg[0] + seg[1]) + 0.5
                thick = float(seg[1] - seg[0] + 1)
            elif series == "actinin_r1p0":
                pred_r = float(y_axis.to_pixel(1.0e-3))
                seg = min(segs, key=lambda s: abs(0.5 * (s[0] + s[1]) - pred_r))
                r_est = 0.5 * (seg[0] + seg[1]) + 0.5
                thick = float(seg[1] - seg[0] + 1)
            else:  # actinin_r2p0
                seg = segs[-1]
                r_est = 0.5 * (seg[0] + seg[1]) + 0.5
                thick = float(seg[1] - seg[0] + 1)
        else:
            pred_r = float(y_axis.to_pixel(10.0 ** (ga + gb * xv)))
            seg = min(segs, key=lambda s: abs(0.5 * (s[0] + s[1]) - pred_r))
            r_est = 0.5 * (seg[0] + seg[1]) + 0.5
            thick = float(seg[1] - seg[0] + 1)
            if abs(r_est - pred_r) > 220:
                out.append({"x": xv, "y": None, "found": False, "note": f"{series} branch occluded at l={xv} um"})
                continue

        half_px = max(0.35 * min(thick, 48.0), 4.0) + y_axis.residual_px
        yv = float(y_axis.to_value(r_est))
        ya, yb = float(y_axis.to_value(r_est - half_px)), float(y_axis.to_value(r_est + half_px))
        out.append(
            {
                "x": xv,
                "y": yv,
                "lo": min(ya, yb),
                "hi": max(ya, yb),
                "px": [c_f, r_est],
                "px_lo": r_est - half_px,
                "px_hi": r_est + half_px,
                "note": f"Fig. S3C {series}; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
            }
        )
    return out
