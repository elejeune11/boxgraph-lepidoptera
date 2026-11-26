"""Shared fast calibration and extraction helpers for liu2008 figure readers.

All functions operate strictly on the RGB numpy array passed into read(image, at);
no PDF rendering or hardcoded pixel coordinates are used.
"""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np
from boxgraph.digitize import Axis


def find_runs(mask: np.ndarray, min_len: int = 1, max_len: int | None = None) -> list[tuple[int, int, float]]:
    """Return contiguous True runs as (start, end_excl, center) in continuous pixel coords."""
    m = np.asarray(mask, dtype=bool)
    if m.size == 0:
        return []
    d = np.diff(np.concatenate(([False], m, [False])).astype(np.int8))
    starts = np.nonzero(d == 1)[0]
    ends = np.nonzero(d == -1)[0]
    out: list[tuple[int, int, float]] = []
    for a, b in zip(starts, ends):
        w = int(b - a)
        if w >= min_len and (max_len is None or w <= max_len):
            out.append((int(a), int(b), 0.5 * (a + b)))
    return out


def dark_neutral_mask(img: np.ndarray, dark_thr: float = 120.0, sat_thr: float = 35.0) -> np.ndarray:
    """Boolean mask of dark neutral (black/dark grey) ink pixels."""
    rgb = img.astype(np.float32)
    gray = rgb.mean(axis=2)
    sat = rgb.max(axis=2) - rgb.min(axis=2)
    return (gray < dark_thr) & (sat < sat_thr)


def pure_black_mask(img: np.ndarray, max_ch: float = 85.0) -> np.ndarray:
    """Boolean mask of black ink strictly darker than red/blue/green fills on all channels."""
    return img.max(axis=2) < max_ch


def find_plot_box(
    img: np.ndarray,
    dark_thr: float = 120.0,
    sat_thr: float = 35.0,
    min_h_frac: float = 0.45,
    min_v_frac: float = 0.45,
) -> dict[str, tuple[int, int, float]]:
    """Detect the closed rectangular plot frame (top, bot, left, right spines).

    Chooses the outermost pair of horizontal lines spanning at least `min_h_frac` of
    image width and vertical lines spanning at least `min_v_frac` of the height
    between top and bottom spines, so text labels outside the box are ignored.
    """
    H, W = img.shape[:2]
    dark = dark_neutral_mask(img, dark_thr=dark_thr, sat_thr=sat_thr)

    # First find long horizontal lines spanning at least min_h_frac * W
    row_cov = dark.mean(axis=1)
    h_runs = find_runs(row_cov >= min_h_frac, min_len=2, max_len=60)
    if len(h_runs) < 2:
        raise RuntimeError(f"Expected >= 2 horizontal box spines, found {h_runs}")
    top, bot = h_runs[0], h_runs[-1]

    # Restrict vertical search to rows between top and bot spines so axis labels don't trigger
    r0, r1 = top[0], bot[1]
    sub = dark[r0:r1, :]
    col_cov = sub.mean(axis=0)
    v_runs = find_runs(col_cov >= max(0.85, min_v_frac), min_len=2, max_len=60)
    if len(v_runs) < 2:
        v_runs = find_runs(col_cov >= min_v_frac, min_len=2, max_len=60)
    if len(v_runs) < 2:
        raise RuntimeError(f"Expected >= 2 vertical box spines in rows [{r0},{r1}], found {v_runs}")
    left, right = v_runs[0], v_runs[-1]
    return {"top": top, "bot": bot, "left": left, "right": right}


def scan_ticks(
    dark: np.ndarray,
    axis: str,
    r0: int,
    r1: int,
    c0: int,
    c1: int,
    min_frac: float = 0.55,
    min_len: int = 2,
    max_len: int = 50,
) -> list[float]:
    """Scan a thin strip parallel to a spine and return tick center coordinates."""
    band = dark[r0:r1, c0:c1]
    if band.size == 0:
        return []
    if axis == "x":
        prof = band.mean(axis=0) >= min_frac
        return [c0 + ctr for _, _, ctr in find_runs(prof, min_len=min_len, max_len=max_len)]
    else:
        prof = band.mean(axis=1) >= min_frac
        return [r0 + ctr for _, _, ctr in find_runs(prof, min_len=min_len, max_len=max_len)]


def select_regular_ticks(cands: Iterable[float], n_expected: int, tol: float = 5.0) -> list[float]:
    """Pick the most evenly spaced subset of length `n_expected` from candidate tick coordinates."""
    c = sorted(float(v) for v in cands)
    if len(c) < n_expected:
        raise RuntimeError(f"Expected {n_expected} regular ticks, got only {len(c)}: {c}")
    best: list[float] | None = None
    best_err = float("inf")
    n = len(c)
    for i in range(n):
        for j in range(i + 1, n):
            step = c[j] - c[i]
            if step < 20.0:
                continue
            seq = [c[i], c[j]]
            ok = True
            for k in range(2, n_expected):
                # Re-estimate spacing from current sequence
                est_step = (seq[-1] - seq[0]) / (len(seq) - 1)
                target = seq[0] + k * est_step
                nearest = min(c, key=lambda v: abs(v - target))
                if abs(nearest - target) > max(tol, 0.02 * est_step):
                    ok = False
                    break
                seq.append(nearest)
            if not ok:
                continue
            xs = np.arange(n_expected, dtype=float)
            ys = np.array(seq, dtype=float)
            slope, intercept = np.polyfit(xs, ys, 1)
            err = float(np.max(np.abs(ys - (intercept + slope * xs))))
            if err < best_err - 1e-6 or (abs(err - best_err) <= 1e-6 and best is not None and (seq[-1] - seq[0]) > (best[-1] - best[0])):
                best_err = err
                best = seq
    if best is None:
        raise RuntimeError(f"Could not match {n_expected} regular ticks in {c}")
    return best


def make_axis(scale: str, marks: list[tuple[float, float]], max_res_px: float = 3.5) -> Axis:
    """Create a BoxGraph Axis and verify calibration misfit."""
    ax = Axis(scale, marks)
    if ax.residual_px > max_res_px:
        raise RuntimeError(f"Calibration residual {ax.residual_px:.2f} px exceeds {max_res_px} px for marks={marks}")
    return ax


def reading_point(
    x_val: float,
    col: float,
    row: float,
    y_axis: Axis,
    half_px: float,
    note: str = "",
    extra_y_px: float = 0.0,
    bar_hi_row: float | None = None,
    bar_lo_row: float | None = None,
    bar_err_px: float | None = None,
    clamp_min: float | None = None,
    clamp_max: float | None = None,
) -> dict:
    """Construct a code-mode reading dict with honest 95% uncertainty interval."""
    h = max(half_px, 0.5) + max(y_axis.residual_px, 0.5) + max(extra_y_px, 0.0)
    y = y_axis.to_value(row)
    a, b = y_axis.to_value(row - h), y_axis.to_value(row + h)
    lo, hi = min(a, b), max(a, b)
    if clamp_min is not None:
        if y < clamp_min and hi >= clamp_min:
            y = clamp_min
        lo = max(clamp_min, lo)
        hi = max(y, hi)
    if clamp_max is not None:
        if y > clamp_max and lo <= clamp_max:
            y = clamp_max
        hi = min(clamp_max, hi)
        lo = min(y, lo)

    rec: dict = {
        "x": float(x_val),
        "y": float(y),
        "lo": float(lo),
        "hi": float(hi),
        "px": [float(col), float(row)],
        "px_lo": float(row - h),
        "px_hi": float(row + h),
        "note": note,
    }
    if bar_hi_row is not None and bar_lo_row is not None:
        be = (bar_err_px if bar_err_px is not None else 2.0) + max(y_axis.residual_px, 0.5)
        r_top = min(bar_hi_row, bar_lo_row)
        r_bot = max(bar_hi_row, bar_lo_row)
        v_hi = y_axis.to_value(r_top)
        v_lo = y_axis.to_value(r_bot)
        v_hi_out = y_axis.to_value(r_top - be)
        v_lo_out = y_axis.to_value(r_bot + be)
        b_lo, b_hi = min(v_lo, v_hi), max(v_lo, v_hi)
        b_lo_o, b_hi_o = min(v_lo_out, v_hi_out), max(v_lo_out, v_hi_out)
        if clamp_min is not None:
            b_lo = max(clamp_min, b_lo)
            b_lo_o = max(clamp_min, b_lo_o)
        if clamp_max is not None:
            b_hi = min(clamp_max, b_hi)
            b_hi_o = min(clamp_max, b_hi_o)
        b_lo = min(b_lo, y)
        b_hi = max(b_hi, y)
        b_lo_o = min(b_lo_o, b_lo)
        b_hi_o = max(b_hi_o, b_hi)
        rec["bar_lo"] = float(b_lo)
        rec["bar_hi"] = float(b_hi)
        rec["bar_lo_outer"] = float(b_lo_o)
        rec["bar_hi_outer"] = float(b_hi_o)
    return rec


def trace_curve_at(
    mask: np.ndarray,
    x_axis: Axis,
    y_axis: Axis,
    at: list[float],
    r_min: int,
    r_max: int,
    half_col: int = 2,
    search_window_px: int = 0,
    min_run: int = 2,
    extra_y_px: float = 0.0,
    clamp_min: float | None = None,
    clamp_max: float | None = None,
    label: str = "",
) -> list[dict]:
    """Read a single-valued curve mask at given x positions with slope-aware intervals."""
    H, W = mask.shape
    r0 = max(0, int(r_min))
    r1 = min(H, int(r_max))

    def sample_col(c_float: float, r_hint: float | None = None) -> tuple[float, float, float] | None:
        """Return (row_center, half_thickness_px, used_col) or None."""
        c0 = int(round(c_float))
        offsets = [0]
        for d in range(1, search_window_px + 1):
            offsets.extend((-d, d))
        for off in offsets:
            cc = c0 + off
            if cc < 0 or cc >= W:
                continue
            lo_c = max(0, cc - half_col)
            hi_c = min(W, cc + half_col + 1)
            col_any = mask[r0:r1, lo_c:hi_c].mean(axis=1) >= 0.4
            runs = find_runs(col_any, min_len=min_run)
            if not runs:
                col_any = mask[r0:r1, lo_c:hi_c].any(axis=1)
                runs = find_runs(col_any, min_len=min_run)
            if not runs:
                continue
            if r_hint is not None and len(runs) > 1:
                chosen = min(runs, key=lambda t: abs((r0 + t[2]) - r_hint))
            else:
                # Pick thickest vertical run (main curve stroke)
                chosen = max(runs, key=lambda t: t[1] - t[0])
            row_ctr = r0 + chosen[2]
            half_thick = 0.5 * (chosen[1] - chosen[0])
            return row_ctr, half_thick, float(cc) + 0.5
        return None

    # First pass to get raw positions for slope estimation
    raw: list[tuple[float, float, tuple[float, float, float] | None]] = []
    for xv in at:
        cf = x_axis.to_pixel(xv)
        raw.append((xv, cf, sample_col(cf)))

    out: list[dict] = []
    for idx, (xv, cf, hit) in enumerate(raw):
        if hit is None:
            out.append({"x": float(xv), "y": None, "note": f"{label}: no curve pixels found at x={xv:g}"})
            continue
        row_ctr, half_thick, used_c = hit
        # Estimate local slope |dr/dc| from neighbors to propagate x-calibration & offset error into y
        slope_px = 0.0
        for nb in (idx - 1, idx + 1):
            if 0 <= nb < len(raw) and raw[nb][2] is not None:
                dc = abs(raw[nb][1] - cf)
                if dc > 1.0:
                    slope_px = max(slope_px, abs(raw[nb][2][0] - row_ctr) / dc)
        col_err = abs(used_c - cf) + max(x_axis.residual_px, 0.5)
        slope_err = slope_px * col_err
        half_px = math.hypot(half_thick, slope_err)
        note = (
            f"{label}: width={2*half_thick:.1f}px, cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px"
            + (f", col_shift={used_c - cf:+.1f}px" if abs(used_c - cf) > 1.5 else "")
        )
        out.append(
            reading_point(
                x_val=xv,
                col=cf,
                row=row_ctr,
                y_axis=y_axis,
                half_px=half_px,
                extra_y_px=extra_y_px,
                note=note,
                clamp_min=clamp_min,
                clamp_max=clamp_max,
            )
        )
    return out
