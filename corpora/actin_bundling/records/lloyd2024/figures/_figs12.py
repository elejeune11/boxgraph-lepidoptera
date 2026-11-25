"""Shared readers for SI Fig. 12A-B AFM profiles (1800 dpi) and SI Fig. 12C ridge height plot (1200 dpi), PDF p. 36."""

from __future__ import annotations

import numpy as np
from scipy import ndimage
from boxgraph import digitize as D
from _common import run_centers, y_interval


def read_afm_profile(image: np.ndarray, at: list[float], xmax: float,
                     color: str = "blue") -> list[dict]:
    """Calibrate SI Fig. 12A or 12B AFM profile axes (0..xmax um, 0..2 um) and trace a colored profile."""
    gray = image.astype(float).mean(axis=2)
    dark = gray < 215.0
    H, W = dark.shape

    c_frac = dark.mean(axis=0)
    c_frac[:200] = 0.0
    r_frac = dark[:, 250:].mean(axis=1)
    c_lines = run_centers(np.where(c_frac > 0.55)[0], gap=6)
    r_lines = run_centers(np.where(r_frac > 0.70)[0], gap=6)
    if len(c_lines) < 2 or len(r_lines) < 2:
        raise RuntimeError(f"Could not detect AFM profile frame: c={c_lines}, r={r_lines}")
    c0, c1 = c_lines[0], c_lines[-1]
    r0, r1 = r_lines[0], r_lines[-1]

    # Check middle y=1 major tick for 3-point y-axis calibration
    yt_mid_cands = run_centers(
        np.where(gray[int(r0 + 0.35 * (r1 - r0)):int(r0 + 0.65 * (r1 - r0)), int(round(c0)) - 42] < 200.0)[0]
        + int(r0 + 0.35 * (r1 - r0)), gap=4
    )
    if len(yt_mid_cands) == 1:
        ya = D.Axis("linear", [(r0, 2.0), (yt_mid_cands[0], 1.0), (r1, 0.0)])
    else:
        ya = D.Axis("linear", [(r0, 2.0), (r1, 0.0)])

    # Check major x-ticks (0, 1, ..., xmax) at offset 45 px below bottom axis
    xt_major = run_centers(
        np.where(gray[int(round(r1)) + 45, max(0, int(c0) - 15):min(W, int(c1) + 15)] < 200.0)[0]
        + max(0, int(c0) - 15), gap=5
    )
    # Match detected major ticks to integer values 0..int(xmax)
    n_int = int(round(xmax))
    step_est = (c1 - c0) / xmax
    x_marks = [(c0, 0.0), (c1, xmax)]
    for k in range(1, n_int):
        pred = c0 + k * step_est
        near = [t for t in xt_major if abs(t - pred) < 20.0]
        if near:
            x_marks.append((min(near, key=lambda t: abs(t - pred)), float(k)))
    x_marks.sort(key=lambda m: m[1])
    xa = D.Axis("linear", x_marks)

    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    in_box = np.zeros((H, W), dtype=bool)
    in_box[int(r0) + 6:int(r1) - 6, int(c0) + 3:int(c1) - 6] = True
    # Exclude top-right legend box
    in_box[int(r0):int(r0 + 0.42 * (r1 - r0)), int(c0 + 0.83 * (c1 - c0)):] = False

    if color == "blue":
        mask = in_box & (B - R > 55.0) & (B - G > 35.0)
    elif color == "green":
        mask = in_box & (G - R > 40.0) & (G - B > 40.0)
    elif color == "magenta":
        mask = in_box & (R - G > 55.0) & (B - G > 35.0)
    else:
        raise ValueError(f"Unsupported profile color {color!r}")

    out = []
    for xv in at:
        col_f = xa.to_pixel(float(xv))
        col = int(round(col_f))
        # Search narrowest window first (col, then +-2 px, then +-5 px) to avoid blending steep slopes
        rows = np.array([], dtype=int)
        for hw in (0, 2, 4):
            c_lo = max(int(c0) + 2, col - hw)
            c_hi = min(int(c1) - 2, col + hw + 1)
            rows = np.where(mask[:, c_lo:c_hi].any(axis=1))[0]
            if len(rows) > 0:
                break
        if len(rows) == 0:
            out.append({"x": float(xv), "y": None, "note": f"no {color} curve pixels at col {col}"})
            continue
        row_c = float(rows.mean() + 0.5)
        half_px = max(0.5 * (rows.max() - rows.min() + 1), 2.5)
        y_val, lo, hi, px_lo, px_hi = y_interval(ya, row_c, half_px)
        out.append({
            "x": float(xv), "y": float(y_val), "lo": float(lo), "hi": float(hi),
            "px": [float(col_f), row_c], "px_lo": px_lo, "px_hi": px_hi,
            "note": f"Profile 4 ({color}), rx={xa.residual_px:.2f} px, ry={ya.residual_px:.2f} px",
        })
    return out


def _calibrate_s12c(image: np.ndarray) -> tuple[D.Axis, float, float, list[float]]:
    gray = image.astype(float).mean(axis=2)
    dark = gray < 100.0
    c0 = run_centers(np.where(dark.mean(axis=0) > 0.50)[0], gap=5)[0]
    r1 = run_centers(np.where(dark.mean(axis=1) > 0.50)[0], gap=5)[-1]
    yt = run_centers(np.where(dark[:int(r1), int(round(c0)) - 20])[0], gap=4)
    xt = run_centers(np.where(dark[int(round(r1)) + 20, int(round(c0)):])[0] + int(round(c0)), gap=4)
    if len(yt) != 3 or len(xt) != 2:
        raise RuntimeError(f"SI Fig. 12C expected 3 y-ticks and 2 x-ticks, got yt={yt}, xt={xt}")
    ya = D.Axis("linear", list(zip(yt, [1400.0, 1200.0, 1000.0])))
    return ya, c0, r1, xt


def read_s12c_mean_sd(image: np.ndarray) -> list[dict]:
    """Read red mean dots and red vertical SD bars for control (x=1) and cyto (x=2) in SI Fig. 12C."""
    ya, c0, r1, xt = _calibrate_s12c(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    red = (R - G > 80.0) & (R - B > 80.0)

    out = []
    labels = ["control", "cyto"]
    for idx, (xc, lbl) in enumerate(zip(xt, labels), start=1):
        c_lo, c_hi = int(round(xc)) - 70, int(round(xc)) + 70
        sub_red = red[:int(r1) - 10, c_lo:c_hi]
        rows_red = np.where(sub_red.any(axis=1))[0]
        widths = sub_red.sum(axis=1)
        stem_w = widths[rows_red].min()
        dot_rows = np.where(widths > stem_w + 15)[0]
        row_mean = float(dot_rows.mean() + 0.5)
        r_top = float(rows_red.min() + 0.5)
        r_bot = float(rows_red.max() + 0.5)

        h_px = 2.0 + max(ya.residual_px, 0.5)
        y_mean, lo, hi, px_lo, px_hi = y_interval(ya, row_mean, 2.0)
        bar_hi = float(ya.to_value(r_top))
        bar_lo = float(ya.to_value(r_bot))
        bar_hi_outer = float(ya.to_value(r_top - h_px))
        bar_lo_outer = float(ya.to_value(r_bot + h_px))
        out.append({
            "x": float(idx),
            "y": float(y_mean), "lo": float(lo), "hi": float(hi),
            "bar_lo": float(min(bar_lo, bar_hi)), "bar_hi": float(max(bar_lo, bar_hi)),
            "bar_lo_outer": float(min(bar_lo_outer, bar_hi_outer)),
            "bar_hi_outer": float(max(bar_lo_outer, bar_hi_outer)),
            "px": [float(xc), row_mean],
            "px_lo": float(r_top), "px_hi": float(r_bot),
            "note": f"{lbl} mean +- SD (axis literally printed um, values are nm; ry={ya.residual_px:.2f} px)",
        })
    return out


def read_s12c_dots(image: np.ndarray, group_index: int) -> list[dict]:
    """Read all 6 individual black dots for group_index (1=control, 2=cyto) in SI Fig. 12C."""
    ya, c0, r1, xt = _calibrate_s12c(image)
    xc = xt[group_index - 1]
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    blk = (R < 40.0) & (G < 40.0) & (B < 40.0)
    c_off = int(round(xc)) - 260
    sub = blk[:int(r1) - 15, c_off:int(round(xc)) + 260]

    lab, n = ndimage.label(sub)
    blobs = []
    for k in range(1, n + 1):
        rr, cc = np.where(lab == k)
        if len(rr) < 1500:
            continue
        w = int(cc.max() - cc.min() + 1)
        h = int(rr.max() - rr.min() + 1)
        blobs.append({
            "r": float(rr.mean() + 0.5),
            "c": float(cc.mean() + 0.5) + c_off,
            "w": w, "h": h, "area": len(rr),
        })

    # Full unsplit dot width is max width (~115-130 px); dots bisected by the vertical red SD line
    # form two half-moons (width < 0.75 * full_w) at identical row centroid. Merge matching pairs.
    full_w = max(b["w"] for b in blobs)
    used = [False] * len(blobs)
    dots = []
    for i, bi in enumerate(blobs):
        if used[i]:
            continue
        if bi["w"] < 0.78 * full_w:
            # Find its twin half-moon
            partner = None
            for j in range(i + 1, len(blobs)):
                if used[j] or blobs[j]["w"] >= 0.78 * full_w:
                    continue
                if abs(blobs[j]["r"] - bi["r"]) < 12.0 and abs(blobs[j]["c"] - bi["c"]) < 1.1 * full_w:
                    partner = j
                    break
            if partner is not None:
                bj = blobs[partner]
                used[i] = used[partner] = True
                tot = bi["area"] + bj["area"]
                rc = (bi["r"] * bi["area"] + bj["r"] * bj["area"]) / tot
                cc = (bi["c"] * bi["area"] + bj["c"] * bj["area"]) / tot
                dots.append((rc, cc))
                continue
        used[i] = True
        dots.append((bi["r"], bi["c"]))

    # Sort ascending by ridge height (descending pixel row)
    dots.sort(key=lambda d: (-d[0], d[1]))
    lbl = "control" if group_index == 1 else "cyto"
    out = []
    for idx, (rc, cc) in enumerate(dots, start=1):
        y_val, lo, hi, px_lo, px_hi = y_interval(ya, rc, 2.5)
        out.append({
            "x": float(idx), "y": float(y_val), "lo": float(lo), "hi": float(hi),
            "px": [float(cc), float(rc)], "px_lo": px_lo, "px_hi": px_hi,
            "note": f"{lbl} dot {idx}/{len(dots)} (ry={ya.residual_px:.2f} px)",
        })
    return out
