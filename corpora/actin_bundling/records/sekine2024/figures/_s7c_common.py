"""Shared calibration and curve-reading helpers for Supplementary Fig. 7c1-c4 (PDF p. 22, 1200 dpi)."""

from __future__ import annotations

import math
import numpy as np
from scipy import ndimage


def _split_runs(idx: np.ndarray, gap: int = 8) -> list[np.ndarray]:
    if len(idx) == 0:
        return []
    cuts = np.where(np.diff(idx) > gap)[0] + 1
    return [g for g in np.split(idx, cuts) if len(g) > 0]


def calibrate_s7c(image: np.ndarray, x_labels: list[float], y_labels_top_to_bottom: list[float]):
    """Detect the rectangular plot box and outward tick centroids on a 1200-dpi Supp. Fig. 7c crop."""
    img = image.astype(float)
    H, W, _ = img.shape
    gray = img.mean(axis=2)
    dark = gray < 120.0

    hr = _split_runs(np.where(dark.sum(axis=1) > W * 0.5)[0], gap=10)
    vc = _split_runs(np.where(dark.sum(axis=0) > H * 0.5)[0], gap=10)
    if len(hr) < 2 or len(vc) < 2:
        raise RuntimeError(f"Could not detect plot box frame (hr={len(hr)}, vc={len(vc)})")

    r_top_in = int(hr[0].max()) + 3
    r_bot_in = int(hr[-1].min()) - 3
    c_left_in = int(vc[0].max()) + 3
    c_right_in = int(vc[-1].min()) - 3

    # Outward x-ticks directly below bottom border
    r0_xt = int(hr[-1].max()) + 6
    r1_xt = min(H, r0_xt + 22)
    x_prof = (255.0 - gray[r0_xt:r1_xt, :]).mean(axis=0)
    x_groups = _split_runs(np.where(x_prof > 100.0)[0], gap=5)
    xt = [float(((g + 0.5) * x_prof[g]).sum() / x_prof[g].sum()) for g in x_groups]
    if len(xt) != len(x_labels):
        raise RuntimeError(f"Expected {len(x_labels)} x-ticks, found {len(xt)}: {xt}")

    # Outward y-ticks directly left of left border
    c1_yt = max(0, int(vc[0].min()) - 6)
    c0_yt = max(0, c1_yt - 22)
    y_prof = (255.0 - gray[:, c0_yt:c1_yt]).mean(axis=1)
    y_groups = _split_runs(np.where(y_prof > 100.0)[0], gap=5)
    yt = [float(((g + 0.5) * y_prof[g]).sum() / y_prof[g].sum()) for g in y_groups]
    if len(yt) != len(y_labels_top_to_bottom):
        raise RuntimeError(f"Expected {len(y_labels_top_to_bottom)} y-ticks, found {len(yt)}: {yt}")

    bx, ax = np.polyfit(xt, x_labels, 1)
    by, ay = np.polyfit(yt, y_labels_top_to_bottom, 1)
    res_x = float(np.max(np.abs(np.array(xt) - (np.array(x_labels) - ax) / bx)))
    res_y = float(np.max(np.abs(np.array(yt) - (np.array(y_labels_top_to_bottom) - ay) / by)))
    if res_x > 2.0 or res_y > 2.0:
        raise RuntimeError(f"Axis calibration residual too high: rx={res_x:.2f}px, ry={res_y:.2f}px")

    return {
        "box": (c_left_in, r_top_in, c_right_in, r_bot_in),
        "ax": float(ax),
        "bx": float(bx),
        "ay": float(ay),
        "by": float(by),
        "res_x": res_x,
        "res_y": res_y,
    }


def read_s7c_fit(image: np.ndarray, at: list[float], x_labels: list[float], y_labels_top_to_bottom: list[float]):
    """Read the dotted black Gaussian fit curve in Supplementary Fig. 7c at requested spatial frequencies."""
    cal = calibrate_s7c(image, x_labels, y_labels_top_to_bottom)
    c0, r0, c1, r1 = cal["box"]
    ax, bx, ay, by = cal["ax"], cal["bx"], cal["ay"], cal["by"]

    sub = image[r0:r1, c0:c1].astype(float)
    R, G, B = sub[:, :, 0], sub[:, :, 1], sub[:, :, 2]
    # Pure neutral dark pixels (dots of the dotted black Gaussian curve, cleanly separated from blue S(q))
    is_black = (R < 100.0) & (G < 100.0) & (B < 100.0) & (np.abs(B - R) < 25.0)
    lab, n = ndimage.label(is_black)

    dots = []
    for k in range(1, n + 1):
        rr, cc = np.nonzero(lab == k)
        if len(rr) < 20:
            continue
        r_c = float(rr.mean()) + 0.5 + r0
        c_c = float(cc.mean()) + 0.5 + c0
        h_half = 0.5 * float(rr.max() - rr.min() + 1)
        x_d = ax + bx * c_c
        y_d = ay + by * r_c
        dots.append((x_d, y_d, c_c, r_c, h_half))
    dots.sort(key=lambda d: d[0])
    if len(dots) < 5:
        raise RuntimeError(f"Too few dotted-curve segments detected: {len(dots)}")

    xs_d = np.array([d[0] for d in dots])
    ys_d = np.array([d[1] for d in dots])
    rs_d = np.array([d[3] for d in dots])
    med_half_h = float(np.median([d[4] for d in dots]))

    out = []
    for xv in at:
        xv = float(xv)
        col = (xv - ax) / bx
        if xv < xs_d[0] - 0.05 or xv > xs_d[-1] + 0.05:
            out.append({"x": xv, "y": None, "note": "outside dotted Gaussian curve range"})
            continue
        # Local quadratic fit to the 5 nearest dot centroids (handles inter-dot gap and peak curvature)
        order = np.argsort(np.abs(xs_d - xv))[:5]
        order = np.sort(order)
        lx, ly, lr = xs_d[order], ys_d[order], rs_d[order]
        coef = np.polyfit(lx - xv, lr, 2)
        r_est = float(coef[2])
        fit_resid_px = float(np.max(np.abs(np.polyval(coef, lx - xv) - lr)))
        y_est = ay + by * r_est

        # Honest 95% pixel half-band: dot half-thickness (~7 px) + calibration misfit + local fit residual + slope*rx
        slope_r_per_c = float(coef[1]) * bx
        err_px = math.hypot(med_half_h * 0.75 + 1.5, abs(slope_r_per_c) * (cal["res_x"] + 1.0)) + cal["res_y"] + fit_resid_px
        err_y = abs(by) * err_px
        out.append({
            "x": xv,
            "y": y_est,
            "lo": max(0.0, y_est - err_y),
            "hi": y_est + err_y,
            "px": [col, r_est],
            "px_lo": r_est - err_px,
            "px_hi": r_est + err_px,
            "note": f"local quadratic over 5 dots (rx={cal['res_x']:.2f}px, ry={cal['res_y']:.2f}px, fit={fit_resid_px:.2f}px)",
        })
    return out


def read_s7c_sim(image: np.ndarray, at: list[float], x_labels: list[float], y_labels_top_to_bottom: list[float]):
    """Read the solid blue simulated S(q) curve and circular bin markers in Supplementary Fig. 7c."""
    cal = calibrate_s7c(image, x_labels, y_labels_top_to_bottom)
    c0, r0, c1, r1 = cal["box"]
    ax, bx, ay, by = cal["ax"], cal["bx"], cal["ay"], cal["by"]

    sub = image[r0:r1, c0:c1].astype(float)
    R, G, B = sub[:, :, 0], sub[:, :, 1], sub[:, :, 2]
    is_blue = (B - R > 30.0) & (B > 70.0)
    # Close tiny internal holes where a black dot sits directly inside a blue marker disk
    blue_closed = ndimage.binary_closing(is_blue, structure=np.ones((9, 9), dtype=bool))
    H_s, W_s = blue_closed.shape

    # Estimate circular marker radius from median vertical span of narrow columns
    col_spans = []
    for xv in at:
        col_s = int(round((float(xv) - ax) / bx)) - c0
        if 0 <= col_s < W_s:
            rr = np.nonzero(blue_closed[:, col_s])[0]
            if len(rr):
                runs = _split_runs(rr, gap=3)
                best = max(runs, key=len)
                col_spans.append(int(best[-1] - best[0] + 1))
    narrow = [s for s in col_spans if s <= 38]
    rad = 0.5 * float(np.median(narrow)) if narrow else 16.5
    dx = min(13, max(8, int(round(0.78 * rad))))
    half_chord = math.sqrt(max(1.0, rad * rad - dx * dx))

    # Extract connected blue chord at col ± dx AND full column span at col for each requested x
    raw = []
    for xv in at:
        xv = float(xv)
        col_f = (xv - ax) / bx
        col_s = int(round(col_f)) - c0
        rr_col = np.nonzero(blue_closed[:, col_s])[0]
        if len(rr_col) == 0:
            raw.append((xv, col_f, None, None, 0, None, None))
            continue
        col_runs = _split_runs(rr_col, gap=4)
        col_best = max(col_runs, key=len)
        c_rmin, c_rmax = int(col_best[0]), int(col_best[-1])

        cL = max(0, col_s - dx)
        cR = min(W_s - 1, col_s + dx)
        both = blue_closed[:, cL] & blue_closed[:, cR] & blue_closed[:, col_s]
        rr = np.nonzero(both)[0]
        if len(rr) == 0:
            rmin, rmax = c_rmin, c_rmax
        else:
            runs = _split_runs(rr, gap=4)
            best = max(runs, key=len)
            rmin, rmax = int(best[0]), int(best[-1])
        raw.append((xv, col_f, rmin, rmax, rmax - rmin + 1, c_rmin, c_rmax))

    mids = [0.5 * (r[5] + r[6]) if r[5] is not None else None for r in raw]
    out = []
    for i, (xv, col_f, rmin, rmax, span, c_rmin, c_rmax) in enumerate(raw):
        if rmin is None:
            out.append({"x": xv, "y": None, "note": "no blue curve found at column"})
            continue
        c_span = c_rmax - c_rmin + 1
        rp = mids[i - 1] if i > 0 and mids[i - 1] is not None else mids[i]
        rn = mids[i + 1] if i + 1 < len(mids) and mids[i + 1] is not None else mids[i]
        is_end = (i == 0 or i == len(raw) - 1)
        if is_end or (span <= 2.0 * half_chord + 5.0 and c_span <= 2.0 * rad + 8.0):
            rc_s = 0.5 * (c_rmin + c_rmax) + 0.5 if is_end else 0.5 * (rmin + rmax) + 0.5
            kind = "mid"
            err_px = 2.0 + cal["res_y"]
        elif rp > c_rmin + rad and rn > c_rmin + rad:
            rc_s = c_rmin + rad
            kind = "peak"
            err_px = 2.5 + cal["res_y"]
        elif rp < c_rmax - rad and rn < c_rmax - rad:
            rc_s = c_rmax - rad
            kind = "trough"
            err_px = 2.5 + cal["res_y"]
        elif span <= 2.0 * half_chord + 5.0:
            rc_s = 0.5 * (rmin + rmax) + 0.5
            kind = "chord-mid"
            err_px = 2.5 + cal["res_y"]
        elif (abs(rp - c_rmin) < abs(rn - c_rmin) and rn > c_rmax) or (abs(rn - c_rmin) < abs(rp - c_rmin) and rp > c_rmax):
            rc_s = c_rmin + rad
            kind = "knee-top"
            err_px = 3.5 + cal["res_y"]
        elif (abs(rp - c_rmax) < abs(rn - c_rmax) and rn < c_rmin) or (abs(rn - c_rmax) < abs(rp - c_rmax) and rp < c_rmin):
            rc_s = c_rmax - rad
            kind = "knee-bot"
            err_px = 3.5 + cal["res_y"]
        else:
            rc_s = 0.5 * (rmin + rmax) + 0.5
            kind = "slope"
            err_px = 0.25 * span + cal["res_y"]

        r_abs = rc_s + r0
        y_est = ay + by * r_abs
        err_y = abs(by) * err_px
        out.append({
            "x": xv,
            "y": y_est,
            "lo": y_est - err_y,
            "hi": y_est + err_y,
            "px": [col_f, r_abs],
            "px_lo": r_abs - err_px,
            "px_hi": r_abs + err_px,
            "note": f"blue marker ({kind}, span={span}px, cspan={c_span}px, rad={rad:.1f}px, rx={cal['res_x']:.2f}px, ry={cal['res_y']:.2f}px)",
        })
    return out
