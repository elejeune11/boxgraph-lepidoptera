"""Shared image calibration and digitization helpers for nastkolb2022 figures (1200 dpi vector OriginLab plots)."""

from __future__ import annotations

import math
from typing import Callable

import numpy as np
from scipy import ndimage


def split_runs(idx: np.ndarray, gap: int = 3) -> list[np.ndarray]:
    if len(idx) == 0:
        return []
    cuts = np.where(np.diff(idx) > gap)[0] + 1
    return [g for g in np.split(idx, cuts) if len(g) > 0]


def weighted_center(indices: np.ndarray, weights: np.ndarray) -> float:
    w = weights[indices]
    s = float(w.sum())
    if s <= 0:
        return float(indices.mean()) + 0.5
    return float(((indices + 0.5) * w).sum() / s)


def find_black_spines(gray: np.ndarray, dark: float = 110.0, min_frac_h: float = 0.35, min_frac_v: float = 0.35):
    """Return (h_spines, v_spines) as lists of (start, end) index tuples for long dark lines."""
    H, W = gray.shape
    d = gray < dark
    h_rows = [(int(r[0]), int(r[-1])) for r in split_runs(np.where(d.sum(axis=1) > min_frac_h * W)[0], gap=6)]
    v_cols = [(int(c[0]), int(c[-1])) for c in split_runs(np.where(d.sum(axis=0) > min_frac_v * H)[0], gap=6)]
    return h_rows, v_cols


def detect_ticks_along_spine(
    mask: np.ndarray,
    spine: tuple[int, int],
    axis: str,
    side: str,
    min_offset: int = 2,
    max_offset: int = 50,
    min_len: int = 10,
    coord_min: int = 0,
    coord_max: int | None = None,
    gap: int = 4,
) -> list[tuple[float, int]]:
    """Detect tick marks protruding from a horizontal ('x') or vertical ('y') spine.

    Returns list of (center_px, length_px) sorted by pixel position.
    """
    H, W = mask.shape
    s0, s1 = spine
    if axis == "x":
        c0 = max(0, coord_min)
        c1 = W if coord_max is None else min(W, coord_max)
        if side == "down":
            r0, r1 = min(H, s1 + min_offset), min(H, s1 + max_offset)
            band = mask[r0:r1, c0:c1]
        else:
            r0, r1 = max(0, s0 - max_offset), max(0, s0 - min_offset)
            band = mask[r0:r1, c0:c1][::-1, :]
        # Contiguous run starting within first 3 pixels from spine
        lens = np.zeros(c1 - c0, dtype=int)
        for j in range(c1 - c0):
            col = band[:, j]
            k = 0
            while k < min(3, len(col)) and not col[k]:
                k += 1
            n = 0
            while k + n < len(col) and col[k + n]:
                n += 1
            lens[j] = n
        runs = split_runs(np.where(lens >= min_len)[0], gap=gap)
        out = []
        for gr in runs:
            if len(gr) <= 28:  # exclude thick box/bar artifacts
                center = float(((gr + c0 + 0.5) * lens[gr]).sum() / lens[gr].sum())
                out.append((center, int(lens[gr].max())))
        return out
    else:
        r0 = max(0, coord_min)
        r1 = H if coord_max is None else min(H, coord_max)
        if side == "left":
            c0, c1 = max(0, s0 - max_offset), max(0, s0 - min_offset)
            band = mask[r0:r1, c0:c1][:, ::-1]
        else:
            c0, c1 = min(W, s1 + min_offset), min(W, s1 + max_offset)
            band = mask[r0:r1, c0:c1]
        lens = np.zeros(r1 - r0, dtype=int)
        for i in range(r1 - r0):
            row = band[i, :]
            k = 0
            while k < min(3, len(row)) and not row[k]:
                k += 1
            n = 0
            while k + n < len(row) and row[k + n]:
                n += 1
            lens[i] = n
        runs = split_runs(np.where(lens >= min_len)[0], gap=gap)
        out = []
        for gr in runs:
            if len(gr) <= 28:
                center = float(((gr + r0 + 0.5) * lens[gr]).sum() / lens[gr].sum())
                out.append((center, int(lens[gr].max())))
        return out


def select_major_ticks(ticks: list[tuple[float, int]], expected: int, frac: float = 0.72) -> list[float]:
    if not ticks:
        raise RuntimeError("No ticks detected")
    longest = max(l for _, l in ticks)
    major = [pos for pos, l in ticks if l >= frac * longest]
    if len(major) != expected:
        raise RuntimeError(f"Expected {expected} major ticks, got {len(major)} (all ticks={ticks})")
    return major


class LinearCal:
    """Linear pixel <-> value calibration: val = a + b * px."""

    def __init__(self, pixels: list[float], values: list[float], max_res_px: float = 2.5):
        if len(pixels) != len(values) or len(pixels) < 2:
            raise RuntimeError(f"Calibration needs >=2 matched ticks, got {len(pixels)} vs {len(values)}")
        p = np.asarray(pixels, dtype=float)
        v = np.asarray(values, dtype=float)
        b, a = np.polyfit(p, v, 1)
        self.a = float(a)
        self.b = float(b)
        self.pixels = p
        self.values = v
        self.res_px = float(np.max(np.abs(p - (v - self.a) / self.b)))
        if self.res_px > max_res_px:
            raise RuntimeError(f"Axis calibration residual {self.res_px:.2f}px exceeds {max_res_px}px (p={pixels}, v={values})")

    def to_val(self, px: float) -> float:
        return self.a + self.b * px

    def to_px(self, val: float) -> float:
        return (val - self.a) / self.b


class BrokenYCal:
    """Two-segment broken linear y-axis (lower segment below break_row, upper segment above break_row)."""

    def __init__(self, lower: LinearCal, upper: LinearCal, break_row: float):
        self.lower = lower
        self.upper = upper
        self.break_row = float(break_row)
        self.res_px = max(lower.res_px, upper.res_px)

    def to_val(self, px: float) -> float:
        return self.upper.to_val(px) if px < self.break_row else self.lower.to_val(px)

    def dval_dpx(self, px: float) -> float:
        return abs(self.upper.b) if px < self.break_row else abs(self.lower.b)


def find_horizontal_ticks_far(
    mask: np.ndarray,
    spine: tuple[int, int],
    side: str,
    far_start: int,
    far_end: int,
    r_min: int,
    r_max: int,
    min_run: int = 8,
) -> list[float]:
    """Find major tick row centers using a column strip far enough from the spine to bypass diagonal axis-break slashes."""
    H, W = mask.shape
    s0, s1 = spine
    if side == "right":
        c0, c1 = min(W, s1 + far_start), min(W, s1 + far_end)
    else:
        c0, c1 = max(0, s0 - far_end), max(0, s0 - far_start)
    band = mask[max(0, r_min) : min(H, r_max), c0:c1]
    counts = band.sum(axis=1)
    runs = split_runs(np.where(counts >= min_run)[0], gap=3)
    out = []
    for gr in runs:
        if len(gr) <= 14:
            out.append(float(((gr + max(0, r_min) + 0.5) * counts[gr]).sum() / counts[gr].sum()))
    return out


def read_origin_boxplots(
    image: np.ndarray,
    x_values: list[float],
    box_color: str,  # "red" or "blue"
    y_cal: LinearCal | BrokenYCal,
    r_plot_top: int,
    r_plot_bot: int,
    c_plot_left: int,
    c_plot_right: int,
    include_whiskers: bool = True,
    min_box_area: int = 800,
    break_gap_max: int = 120,
    x_cal: LinearCal | None = None,
) -> list[dict]:
    """Read OriginLab red or blue boxplots (median, 25-75th box, 1-99th whiskers).

    Works for both unbroken and broken y-axes (bridging the white break gap along the whisker stem).
    """
    im = image.astype(float)
    H, W, _ = im.shape
    R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]
    gray = im.mean(axis=2)

    if box_color == "red":
        cmask = (R > 120) & (R - G > 75) & (R - B > 75)
    elif box_color == "blue":
        cmask = (B > 65) & (B - R > 42) & (G > 35)
    else:
        raise ValueError(f"Unknown box_color {box_color!r}")

    # Restrict strictly inside plot area (excludes axis labels and titles)
    sub_mask = np.zeros((H, W), dtype=bool)
    sub_mask[r_plot_top:r_plot_bot, c_plot_left:c_plot_right] = cmask[r_plot_top:r_plot_bot, c_plot_left:c_plot_right]

    # Close across the black median line splitting the colored box into upper & lower halves!
    closed = ndimage.binary_closing(sub_mask, structure=np.ones((25, 5), dtype=bool))
    lab, n = ndimage.label(closed)

    boxes = []
    for k in range(1, n + 1):
        rr, cc = np.nonzero(lab == k)
        if len(rr) < min_box_area:
            continue
        c0, c1 = int(cc.min()), int(cc.max())
        r0, r1 = int(rr.min()), int(rr.max())
        # Filter out tall skinny noise
        if (c1 - c0) < 18:
            continue
        boxes.append((0.5 * (c0 + c1) + 0.5, c0, c1, r0, r1))

    boxes.sort(key=lambda b: b[0])
    if len(boxes) != len(x_values):
        raise RuntimeError(f"Expected {len(x_values)} {box_color} boxes, found {len(boxes)}: {boxes}")

    dark = np.max(im, axis=2) < 85.0
    out = []
    for idx, (c_mid, c0, c1, r_fill_top, r_fill_bot) in enumerate(boxes):
        bw = c1 - c0 + 1
        # Inner columns excluding vertical box walls
        ic0 = c0 + max(6, int(0.18 * bw))
        ic1 = c1 - max(6, int(0.18 * bw))

        # Locate exact outer black top/bottom border of the box and interior median bar.
        # Inside a pure red/blue box, anti-aliased thin black lines have max(R,G,B) well below fill brightness (<150) and are NOT box fill!
        bar_ink = (np.max(im, axis=2) < 150.0) & (~cmask)
        scan_r0 = max(r_plot_top, r_fill_top - 14)
        scan_r1 = min(r_plot_bot, r_fill_bot + 14)
        horiz_frac = bar_ink[scan_r0:scan_r1, ic0 : ic1 + 1].mean(axis=1)
        h_bars = split_runs(np.where(horiz_frac > 0.65)[0], gap=2)
        bar_rows = [float(gr.mean()) + 0.5 + scan_r0 for gr in h_bars if len(gr) <= 18]

        # Top border touches r_fill_top; bottom border touches r_fill_bot; never grab an external whisker cap!
        top_cands = [r for r in bar_rows if r_fill_top - 12 <= r <= r_fill_top + 4]
        bot_cands = [r for r in bar_rows if r_fill_bot - 4 <= r <= r_fill_bot + 12]
        r_box_top = min(top_cands, key=lambda r: abs(r - (r_fill_top - 4))) if top_cands else float(r_fill_top)
        r_box_bot = min(bot_cands, key=lambda r: abs(r - (r_fill_bot + 4))) if bot_cands else float(r_fill_bot)
        interior_bars = [r for r in bar_rows if r > max(r_box_top + 4.0, r_fill_top) and r < min(r_box_bot - 4.0, r_fill_bot)]
        extra_med_px = 0.0
        if interior_bars:
            # Pick the interior bar with darkest minimum brightness across inner columns
            r_med = min(interior_bars, key=lambda r: float(im[int(round(r)), ic0 : ic1 + 1].max(axis=1).mean()))
        else:
            # Strictly check inside the colored fill [r_fill_top+1, r_fill_bot-1] for a thin anti-aliased median line
            if r_fill_bot - r_fill_top > 6:
                i0, i1 = int(r_fill_top) + 1, int(r_fill_bot)
                prof = im[i0:i1, ic0 : ic1 + 1].max(axis=2).mean(axis=1)
                if len(prof) > 0 and float(prof.min()) < 130.0:
                    r_med = float(i0 + int(np.argmin(prof))) + 0.5
                else:
                    # Single unbroken fill band (median merged into Q25 or Q75 border in narrow box):
                    # place center at box midpoint and widen reading interval to span full IQR [r_box_top, r_box_bot]!
                    r_med = 0.5 * (r_box_top + r_box_bot)
                    extra_med_px = 0.5 * (r_box_bot - r_box_top)
            else:
                r_med = 0.5 * (r_box_top + r_box_bot)
                extra_med_px = 0.5 * (r_box_bot - r_box_top)

        # Trace vertical whisker stem upward from r_box_top - 2 and downward from r_box_bot + 2
        # Check narrow central band around c_mid (stem is ~4-8 px wide)
        ci = int(round(c_mid))
        stem_band = dark[:, max(0, ci - 6) : min(W, ci + 7)].any(axis=1)
        # Also check cap width (horizontal whisker cap is ~40-60% of box width, centered on ci)
        cap_half = max(12, int(0.28 * bw))
        cap_frac = dark[:, max(0, ci - cap_half) : min(W, ci + cap_half + 1)].mean(axis=1)

        # Walk upward from r_box_top - 3, bridging breaks up to break_gap_max if broken axis
        max_up_gap = break_gap_max if isinstance(y_cal, BrokenYCal) else 8
        r_up = int(math.floor(r_box_top)) - 2
        last_dark = int(math.floor(r_box_top))
        gap_cnt = 0
        while r_up >= r_plot_top:
            if stem_band[r_up]:
                last_dark = r_up
                gap_cnt = 0
            else:
                gap_cnt += 1
                if gap_cnt > max_up_gap:
                    break
            r_up -= 1
        # Refine upper whisker cap: look for horizontal cap within 12 px below last_dark
        cap_search_u = np.where(cap_frac[max(r_plot_top, last_dark - 4) : min(r_plot_bot, last_dark + 15)] > 0.65)[0]
        if len(cap_search_u) > 0:
            r_p99 = float(cap_search_u.mean()) + 0.5 + max(r_plot_top, last_dark - 4)
        else:
            r_p99 = float(last_dark) + 0.5

        # Walk downward from r_box_bot + 2
        r_dn = int(math.ceil(r_box_bot)) + 2
        last_dark_dn = int(math.ceil(r_box_bot))
        gap_cnt = 0
        while r_dn < r_plot_bot:
            if stem_band[r_dn]:
                last_dark_dn = r_dn
                gap_cnt = 0
            else:
                gap_cnt += 1
                if gap_cnt > 8:
                    break
            r_dn += 1
        cap_search_d = np.where(cap_frac[max(r_plot_top, last_dark_dn - 14) : min(r_plot_bot, last_dark_dn + 5)] > 0.65)[0]
        if len(cap_search_d) > 0:
            r_p01 = float(cap_search_d.mean()) + 0.5 + max(r_plot_top, last_dark_dn - 14)
        else:
            r_p01 = float(last_dark_dn) + 0.5

        y_med = y_cal.to_val(r_med)
        y_q25 = y_cal.to_val(r_box_bot)
        y_q75 = y_cal.to_val(r_box_top)
        y_p01 = min(y_q25, y_cal.to_val(r_p01))
        y_p99 = max(y_q75, y_cal.to_val(r_p99))

        dydp_med = y_cal.dval_dpx(r_med) if isinstance(y_cal, BrokenYCal) else abs(y_cal.b)
        dydp_hi = y_cal.dval_dpx(r_p99) if isinstance(y_cal, BrokenYCal) else abs(y_cal.b)
        dydp_lo = y_cal.dval_dpx(r_p01) if isinstance(y_cal, BrokenYCal) else abs(y_cal.b)

        err_med_px = 2.5 + y_cal.res_px + extra_med_px
        err_med = dydp_med * err_med_px
        pad_hi = dydp_hi * (3.0 + y_cal.res_px)
        pad_lo = dydp_lo * (3.0 + y_cal.res_px)

        xv = float(x_values[idx])
        rec = {
            "x": xv,
            "y": y_med,
            "lo": y_med - err_med,
            "hi": y_med + err_med,
            "px": [c_mid, r_med],
            "note": (
                f"boxplot median={y_med:.4g}, IQR=[{y_q25:.4g}, {y_q75:.4g}], "
                f"1-99th=[{y_p01:.4g}, {y_p99:.4g}], ry={y_cal.res_px:.2f}px"
            ),
        }
        if include_whiskers:
            rec.update({
                "bar_lo": y_p01,
                "bar_hi": y_p99,
                "bar_lo_outer": y_p01 - pad_lo,
                "bar_hi_outer": y_p99 + pad_hi,
                "px_lo": r_p99,
                "px_hi": r_p01,
            })
        else:
            rec.update({
                "px_lo": r_med - err_med_px,
                "px_hi": r_med + err_med_px,
            })
        out.append(rec)
    return out


def read_colored_curve_at(
    mask: np.ndarray,
    x_cal: LinearCal,
    y_cal: LinearCal,
    at: list[float],
    r_top: int,
    r_bot: int,
    c_left: int,
    c_right: int,
    search_cols: int = 8,
    window_half: int = 2,
    y_hint_fn: Callable[[float], tuple[float, float]] | None = None,
) -> list[dict]:
    """Read a continuous curve given its clean color mask at requested x values."""
    H, W = mask.shape
    sub = np.zeros_like(mask)
    sub[r_top:r_bot, c_left:c_right] = mask[r_top:r_bot, c_left:c_right]

    out = []
    for xv in at:
        pc = x_cal.to_px(xv)
        col = int(round(pc))
        # Restrict row range if hint provided (to avoid legend box or distant crossing artifacts)
        rl, rh = r_top, r_bot
        if y_hint_fn is not None:
            ylo, yhi = y_hint_fn(xv)
            r1 = y_cal.to_px(yhi)
            r2 = y_cal.to_px(ylo)
            rl = max(r_top, int(math.floor(min(r1, r2))))
            rh = min(r_bot, int(math.ceil(max(r1, r2))))

        # Collect runs in small column neighborhood around col
        found_row = None
        half_h = 2.0
        used_col = col
        for dc in sorted(range(-search_cols, search_cols + 1), key=abs):
            c = col + dc
            if not (c_left <= c < c_right):
                continue
            rows = np.nonzero(sub[rl:rh, c])[0] + rl
            if len(rows) == 0:
                continue
            grs = split_runs(rows, gap=4)
            # Pick the most compact run if multiple (or median if close)
            if len(grs) == 1:
                gr = grs[0]
            else:
                # Prefer run closest to neighbors in c ± 10
                neigh_rows = []
                for k in (-12, -8, -4, 4, 8, 12):
                    ck = c + k
                    if c_left <= ck < c_right:
                        nr = np.nonzero(sub[rl:rh, ck])[0] + rl
                        if len(nr) > 0:
                            neigh_rows.append(float(np.median(nr)))
                ref = float(np.median(neigh_rows)) if neigh_rows else float(np.median(rows))
                gr = min(grs, key=lambda g: abs((g.mean() + 0.5) - ref))
            found_row = float(gr.mean()) + 0.5
            half_h = 0.5 * (float(gr[-1] - gr[0] + 1))
            used_col = c
            break

        if found_row is None:
            out.append({"x": float(xv), "y": None, "found": False, "note": f"no curve pixels near x={xv:g}"})
            continue

        # Refine with median across [used_col - window_half, used_col + window_half] near found_row to suppress pixel noise
        local_centers = []
        local_cols = []
        for c in range(max(c_left, used_col - window_half), min(c_right, used_col + window_half + 1)):
            rows = np.nonzero(sub[max(rl, int(found_row - 40)) : min(rh, int(found_row + 40)), c])[0] + max(
                rl, int(found_row - 40)
            )
            if len(rows) > 0:
                local_centers.append(float(rows.mean()) + 0.5)
                local_cols.append(float(c) + 0.5)
        if len(local_centers) >= 2:
            # Linear fit across local window to evaluate at exact subpixel pc and get local slope
            slope, intercept = np.polyfit(local_cols, local_centers, 1)
            r_eval = float(intercept + slope * pc) if abs(pc - (used_col + 0.5)) <= search_cols else float(np.median(local_centers))
        else:
            slope = 0.0
            r_eval = found_row

        # Uncertainty: line half-thickness + 1 px + calibration residuals + slope*(x_offset + x_res)
        offset_px = abs(pc - (used_col + 0.5))
        err_px = half_h + 1.0 + y_cal.res_px + abs(slope) * (min(offset_px, 2.0) + x_cal.res_px)
        y_val = y_cal.to_val(r_eval)
        ya = y_cal.to_val(r_eval - err_px)
        yb = y_cal.to_val(r_eval + err_px)
        out.append({
            "x": float(xv),
            "y": y_val,
            "lo": min(ya, yb),
            "hi": max(ya, yb),
            "px": [pc, r_eval],
            "px_lo": r_eval - err_px,
            "px_hi": r_eval + err_px,
            "note": f"rx={x_cal.res_px:.2f}px, ry={y_cal.res_px:.2f}px, half_h={half_h:.1f}px",
        })
    return out
