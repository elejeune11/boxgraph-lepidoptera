"""Shared axis-calibration and curve/marker digitization helpers for shin2009.

All routines operate directly on the rendered RGB crop array passed into
`read(image, at)` — zero PDF I/O or external rendering inside reader code.
"""
from __future__ import annotations

import math
import numpy as np
from boxgraph.digitize import Axis


def runs_1d(idx: np.ndarray) -> list[tuple[int, int, float]]:
    """Return contiguous (start, end, center) runs in continuous coordinates."""
    if len(idx) == 0:
        return []
    out: list[tuple[int, int, float]] = []
    s = p = int(idx[0])
    for v in idx[1:]:
        v = int(v)
        if v == p + 1:
            p = v
        else:
            out.append((s, p, 0.5 * (s + p) + 0.5))
            s = p = v
    out.append((s, p, 0.5 * (s + p) + 0.5))
    return out


def value_interval(axis: Axis, px: float, half_px: float) -> tuple[float, float, float]:
    """Return (y, lo, hi) in axis units for pixel coordinate `px +/- half_px`."""
    y = float(axis.to_value(px))
    a = float(axis.to_value(px - half_px))
    b = float(axis.to_value(px + half_px))
    return y, min(a, b), max(a, b)


def find_box_frame(image: np.ndarray, frac: float = 0.45, thresh: float = 100.0):
    """Locate outer box frame lines (top, bottom, left, right) as (start, end, center)."""
    gray = image.astype(float).mean(axis=2)
    dark = gray < thresh
    H, W = dark.shape
    h_lines = runs_1d(np.where(dark.sum(axis=1) > frac * W)[0])
    v_lines = runs_1d(np.where(dark.sum(axis=0) > frac * H)[0])
    if len(h_lines) < 2 or len(v_lines) < 2:
        raise RuntimeError(f"Expected closed plot frame, found H={len(h_lines)} V={len(v_lines)}")
    return h_lines[0], h_lines[-1], v_lines[0], v_lines[-1]


def black_mask(image: np.ndarray, gray_max: float = 90.0, sat_max: float = 35.0) -> np.ndarray:
    rgb = image.astype(float)
    gray = rgb.mean(axis=2)
    sat = rgb.max(axis=2) - rgb.min(axis=2)
    return (gray < gray_max) & (sat < sat_max)


# ---------------------------------------------------------------------------
# Panel calibrations
# ---------------------------------------------------------------------------

def calibrate_fig1c(image: np.ndarray):
    """Calibrate Fig. 1(c) (fascin:actin I(q) vs q) from matched frame borders and inward ticks."""
    top, bot, left, right = find_box_frame(image)
    rT, rB = top[1], bot[0]
    cL, cR = left[1], right[0]
    blk = black_mask(image)

    bt = [m for _, _, m in runs_1d(np.where(blk[rB - 25 : rB - 12, :].mean(axis=0) > 0.5)[0])
          if cL + 20 < m < cR - 20]
    tt = [m for _, _, m in runs_1d(np.where(blk[rT + 12 : rT + 25, :].mean(axis=0) > 0.5)[0])
          if cL + 20 < m < cR - 20]
    x_matched = [c for c in bt if any(abs(c - ct) < 4.0 for ct in tt)]
    if len(x_matched) != 6:
        raise RuntimeError(f"Fig. 1(c) expected 6 matched x-ticks (0.11..0.135), got {len(x_matched)}: {x_matched}")
    x_marks = list(zip(x_matched, [0.11, 0.115, 0.12, 0.125, 0.13, 0.135])) + [(right[2], 0.14)]
    x_axis = Axis("linear", x_marks)

    lt = [m for _, _, m in runs_1d(np.where(blk[:, cL + 12 : cL + 25].mean(axis=1) > 0.5)[0])
          if rT + 20 < m < rB - 20]
    rt = [m for _, _, m in runs_1d(np.where(blk[:, cR - 25 : cR - 12].mean(axis=1) > 0.5)[0])
          if rT + 20 < m < rB - 20]
    y_matched = [0.5 * (yl + min(rt, key=lambda yr: abs(yr - yl)))
                 for yl in lt if any(abs(yl - yr) < 4.0 for yr in rt)]
    if len(y_matched) != 2:
        raise RuntimeError(f"Fig. 1(c) expected 2 matched y-ticks (2, 1), got {len(y_matched)}: {y_matched}")
    y_axis = Axis("linear", [(top[2], 3.0), (y_matched[0], 2.0), (y_matched[1], 1.0), (bot[2], 0.0)])
    return x_axis, y_axis, (cL, rT, cR, rB)


def calibrate_fig1d(image: np.ndarray, region: str):
    """Calibrate Fig. 1(d) (espin:actin I(q) vs q) for 'upper' ([100,200]) or 'lower' ([40,45]) panel."""
    top, bot, left, right = find_box_frame(image)
    rT, rB = top[1], bot[0]
    cL, cR = left[1], right[0]
    blk = black_mask(image)

    # Top and bottom x-ticks come from the upper and lower sub-panels joined at the '//' break.
    if region == "upper":
        xt = [m for _, _, m in runs_1d(np.where(blk[rT + 12 : rT + 25, :].mean(axis=0) > 0.5)[0])
              if cL + 20 < m < cR - 20]
    elif region == "lower":
        xt = [m for _, _, m in runs_1d(np.where(blk[rB - 25 : rB - 12, :].mean(axis=0) > 0.5)[0])
              if cL + 20 < m < cR - 20]
    else:
        raise ValueError(f"Unknown Fig. 1(d) region {region!r}")
    if len(xt) != 6:
        raise RuntimeError(f"Fig. 1(d) {region} expected 6 x-ticks (0.11..0.135), got {len(xt)}: {xt}")
    x_axis = Axis("linear", list(zip(xt, [0.11, 0.115, 0.12, 0.125, 0.13, 0.135])))

    # Detect left and right inward black ticks; select the two long major ticks (100 and 45)
    # via tick length >= 50 px (excludes the shorter ~39 px '//' break ticks and curve crossings).
    def long_ticks(side: str) -> list[float]:
        out = []
        sl = blk[:, cL + 5 : cL + 20] if side == "L" else blk[:, cR - 20 : cR - 5]
        for s, e, m in runs_1d(np.where(sl.mean(axis=1) > 0.5)[0]):
            if not (rT + 20 < m < rB - 20):
                continue
            w = 0
            if side == "L":
                while cL + 1 + w < cR and blk[s : e + 1, cL + 1 + w].mean() > 0.4:
                    w += 1
            else:
                while cR - 1 - w > cL and blk[s : e + 1, cR - 1 - w].mean() > 0.4:
                    w += 1
            # True major y-ticks have stroke height ~28-36 px and inward length >= 50 px
            if 26 <= (e - s + 1) <= 40 and w >= 50:
                out.append(m)
        return out

    lt = long_ticks("L")
    rt = long_ticks("R")
    if len(lt) != 2 or len(rt) != 2:
        raise RuntimeError(f"Fig. 1(d) expected 2 long left/right y-ticks (100, 45), got L={lt} R={rt}")
    lr_diff_px = max(abs(lt[0] - rt[0]), abs(lt[1] - rt[1]))
    # Use left tick positions (adjacent to printed tick labels) and propagate left-right offset into error
    if region == "upper":
        y_axis = Axis("linear", [(top[2], 200.0), (lt[0], 100.0)])
    else:
        y_axis = Axis("linear", [(lt[1], 45.0), (bot[2], 40.0)])
    return x_axis, y_axis, (cL, rT, cR, rB), lr_diff_px


def calibrate_fig1e(image: np.ndarray):
    """Calibrate Fig. 1(e) (overtwist [deg] vs R) from frame borders and inward ticks."""
    top, bot, left, right = find_box_frame(image)
    rT, rB = top[1], bot[0]
    cL, cR = left[1], right[0]
    blk = black_mask(image)

    # Bottom x-ticks (R = 0.10, 0.15, 0.20, 0.25, 0.30, 0.35; R = 0.05 is occluded by red circle at y=0,
    # and the right frame border sits outside R = 0.40 by half a stroke width)
    bt = [m for s, e, m in runs_1d(np.where(blk[rB - 25 : rB - 12, :].mean(axis=0) > 0.5)[0])
          if cL + 350 < m < cR - 20 and (e - s + 1) >= 22]
    if len(bt) != 6:
        raise RuntimeError(f"Fig. 1(e) expected 6 unobstructed bottom x-ticks (0.1..0.35), got {len(bt)}: {bt}")
    x_marks = [(left[2], 0.0)] + list(zip(bt, [0.10, 0.15, 0.20, 0.25, 0.30, 0.35]))
    x_axis = Axis("linear", x_marks)

    # Left y-ticks (1.0, 0.8, 0.6, 0.4, 0.2) plus top (1.2) and bottom (0.0)
    lt = [m for s, e, m in runs_1d(np.where(blk[:, cL + 12 : cL + 25].mean(axis=1) > 0.5)[0])
          if rT + 20 < m < rB - 20 and (e - s + 1) >= 20]
    if len(lt) != 5:
        raise RuntimeError(f"Fig. 1(e) expected 5 interior left y-ticks (1.0..0.2), got {len(lt)}: {lt}")
    y_marks = [(top[2], 1.2)] + list(zip(lt, [1.0, 0.8, 0.6, 0.4, 0.2])) + [(bot[2], 0.0)]
    y_axis = Axis("linear", y_marks)
    return x_axis, y_axis, (cL, rT, cR, rB)


def calibrate_fig3b(image: np.ndarray):
    """Calibrate Fig. 3(b) (overtwist [rad] vs (beta*U*rho)^2) strictly from the 5 major x-ticks
    (0.05..0.25) and 6 major y-ticks (0.02..0.12), avoiding the hand-placed black x-axis arrow offset."""
    gray = image.astype(float).mean(axis=2)
    dark = gray < 100.0
    H, W = dark.shape
    h_lines = runs_1d(np.where(dark.sum(axis=1) > 0.45 * W)[0])
    v_lines = runs_1d(np.where(dark.sum(axis=0) > 0.45 * H)[0])
    if not h_lines or not v_lines:
        raise RuntimeError("Fig. 3(b): could not locate main axes")
    r_ax = h_lines[0][0]
    c_ax = v_lines[0][1]

    # Major x-ticks stick >= 22 px above r_ax (minor ticks ~13 px)
    x_maj = []
    for s, e, m in runs_1d(np.where(dark[r_ax - 10 : r_ax - 2, :].mean(axis=0) > 0.5)[0]):
        if not (c_ax + 100 < m < W - 300):
            continue
        h = 0
        while r_ax - 1 - h > 0 and dark[r_ax - 1 - h, s : e + 1].mean() > 0.4:
            h += 1
        if h >= 22:
            x_maj.append(m)
    if len(x_maj) != 5:
        raise RuntimeError(f"Fig. 3(b) expected 5 major x-ticks (0.05..0.25), got {len(x_maj)}: {x_maj}")
    x_axis = Axis("linear", list(zip(x_maj, [0.05, 0.10, 0.15, 0.20, 0.25])))

    # Major y-ticks stick 25..40 px right of c_ax with stroke height 15..28 px (excludes arrowhead)
    y_maj = []
    for s, e, m in runs_1d(np.where(dark[:, c_ax + 2 : c_ax + 10].mean(axis=1) > 0.5)[0]):
        if not (120 < m < r_ax - 50):
            continue
        if not (15 <= (e - s + 1) <= 28):
            continue
        w = 0
        while c_ax + 1 + w < W and dark[s : e + 1, c_ax + 1 + w].mean() > 0.4:
            w += 1
        if 25 <= w <= 42:
            y_maj.append(m)
    if len(y_maj) != 6:
        raise RuntimeError(f"Fig. 3(b) expected 6 major y-ticks (0.12..0.02), got {len(y_maj)}: {y_maj}")
    y_axis = Axis("linear", list(zip(y_maj, [0.12, 0.10, 0.08, 0.06, 0.04, 0.02])))
    return x_axis, y_axis, (c_ax, r_ax)


# ---------------------------------------------------------------------------
# Generic curve reader helper
# ---------------------------------------------------------------------------

def sample_curve_from_mask(
    mask: np.ndarray,
    x_axis: Axis,
    y_axis: Axis,
    at: list[float],
    r_lo: int,
    r_hi: int,
    extra_err_px: float = 0.0,
    search_dx: int = 12,
    note_prefix: str = "",
) -> list[dict]:
    """Sample a single-valued curve from a boolean mask at x positions `at`.

    If the exact column `c` is briefly occluded (e.g., by a vertical dashed guide line),
    interpolates cleanly from the nearest clean mask columns within +/- `search_dx`.
    """
    H, W = mask.shape
    cal_y = max(y_axis.residual_px, 0.5)
    out: list[dict] = []

    def col_span(col: int) -> tuple[float, float] | None:
        if not (0 <= col < W):
            return None
        rows = np.where(mask[r_lo:r_hi, col])[0] + r_lo
        if len(rows) < 3:
            return None
        # Pick the thickest contiguous vertical run (rejects isolated 1-2 px speckles)
        rs = runs_1d(rows)
        s, e, ctr = max(rs, key=lambda t: t[1] - t[0])
        return ctr, float(e - s + 1)

    for xv in at:
        c_f = float(x_axis.to_pixel(xv))
        c0 = int(round(c_f))
        hit = col_span(c0)
        interp_note = ""
        if hit is not None:
            r_c, thick = hit
            # Also check +/- 3 px for local slope contribution
            h_l = col_span(c0 - 3)
            h_r = col_span(c0 + 3)
            slope_px = abs(h_r[0] - h_l[0]) / 6.0 if (h_l and h_r) else 0.0
            half_px = 0.5 * thick + 1.0 + cal_y + extra_err_px + 0.5 * slope_px
        else:
            # Interpolate across brief occlusion (such as vertical dashed line)
            left_h = next(((c0 - d, col_span(c0 - d)) for d in range(1, search_dx + 1) if col_span(c0 - d)), None)
            right_h = next(((c0 + d, col_span(c0 + d)) for d in range(1, search_dx + 1) if col_span(c0 + d)), None)
            if left_h is None or right_h is None:
                out.append({"x": float(xv), "y": None, "note": f"{note_prefix}no curve pixels near col {c0}"})
                continue
            cl, (rl, tl) = left_h
            cr, (rr, tr) = right_h
            t = (c_f - cl) / (cr - cl)
            r_c = rl + t * (rr - rl)
            thick = 0.5 * (tl + tr)
            half_px = 0.5 * thick + 2.0 + cal_y + extra_err_px + 0.25 * abs(rr - rl)
            interp_note = f"; bridged gap [{cl},{cr}]"

        y_val, y_lo, y_hi = value_interval(y_axis, r_c, half_px)
        out.append({
            "x": float(xv),
            "y": y_val,
            "lo": y_lo,
            "hi": y_hi,
            "px": [c_f, r_c],
            "px_lo": r_c - half_px,
            "px_hi": r_c + half_px,
            "note": f"{note_prefix}thick={thick:.0f}px, cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px{interp_note}",
        })
    return out


def track_black_curve(
    blk: np.ndarray,
    c_start: int,
    c_end: int,
    r_seed: int,
    r_min: int,
    r_max: int,
    exclude_cols: list[tuple[int, int]] | None = None,
    max_step: int = 30,
    max_thick: int = 65,
) -> dict[int, tuple[float, float]]:
    """Track a single continuous black curve from `c_start` to `c_end` starting at row `r_seed`.

    Returns dict col -> (center_row, thickness_px), cleanly skipping excluded columns
    (e.g., vertical dashed lines) and rejecting black text/arrows outside the tracked band.
    """
    excl = exclude_cols or []

    def is_excl(c: int) -> bool:
        return any(a <= c <= b for a, b in excl)

    tracked: dict[int, tuple[float, float]] = {}
    cur = float(r_seed)
    gap = 0
    for c in range(c_start, c_end + 1):
        if is_excl(c):
            gap += 1
            continue
        win = min(max_step + gap * 4, 140)
        lo = max(r_min, int(cur - win))
        hi = min(r_max, int(cur + win))
        rows = np.where(blk[lo:hi, c])[0] + lo
        if len(rows) < 4:
            gap += 1
            continue
        rs = runs_1d(rows)
        # Filter runs that look like a single curve stroke (not an arrow or vertical shaft)
        cands = [(s, e, ctr) for s, e, ctr in rs if 6 <= (e - s + 1) <= max_thick]
        if not cands:
            gap += 1
            continue
        s, e, ctr = min(cands, key=lambda t: abs(t[2] - cur))
        if abs(ctr - cur) > win:
            gap += 1
            continue
        tracked[c] = (ctr, float(e - s + 1))
        cur = ctr
        gap = 0
    return tracked
