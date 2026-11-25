"""Shared image-calibration and curve/marker digitization helpers for barberi2023.

All functions operate directly on the RGB numpy array `img` passed to `read(img, at)`
by BoxGraph (no PDF rendering, no external I/O, no hardcoded pixel coordinates).
"""

from __future__ import annotations

import math
import numpy as np
from scipy import ndimage


def _groups(indices: np.ndarray, max_gap: int = 3) -> list[tuple[int, int, float]]:
    """Group sorted 1D integer indices into contiguous runs (start, end, center)."""
    if len(indices) == 0:
        return []
    splits = np.where(np.diff(indices) > max_gap)[0] + 1
    out = []
    for g in np.split(indices, splits):
        out.append((int(g[0]), int(g[-1]), float(g.mean()) + 0.5))
    return out


def pink_mask(img: np.ndarray) -> np.ndarray:
    """Magenta/pink ink used for velocity V, H reference, and advective flux curves."""
    R = img[:, :, 0].astype(np.int16)
    G = img[:, :, 1].astype(np.int16)
    B = img[:, :, 2].astype(np.int16)
    return (R - G > 50) & (R > 130) & (B > 50)


def dark_mask(img: np.ndarray, thresh: float = 135.0) -> np.ndarray:
    """Dark neutral (black/charcoal) ink, excluding colored pink/green/blue curves."""
    R = img[:, :, 0].astype(np.float32)
    G = img[:, :, 1].astype(np.float32)
    B = img[:, :, 2].astype(np.float32)
    gray = (R + G + B) / 3.0
    spread = np.maximum(np.maximum(R, G), B) - np.minimum(np.minimum(R, G), B)
    return (gray < thresh) & (spread < 35.0)


def gray_curve_mask(img: np.ndarray, lo: float = 75.0, hi: float = 180.0) -> np.ndarray:
    """Medium-gray neutral ink used for diffusive flux curves in SM Fig. 5(b)-(d)."""
    R = img[:, :, 0].astype(np.float32)
    G = img[:, :, 1].astype(np.float32)
    B = img[:, :, 2].astype(np.float32)
    gray = (R + G + B) / 3.0
    spread = np.maximum(np.maximum(R, G), B) - np.minimum(np.minimum(R, G), B)
    return (gray >= lo) & (gray <= hi) & (spread < 20.0)


def find_frame(img: np.ndarray, include_pink_right: bool = False, min_frac: float = 0.50) -> dict:
    """Locate the 4 rectangular plot spines (c0, r0, c1, r1) and spine spans."""
    dark = dark_mask(img, thresh=140.0)
    spine_mask = dark | pink_mask(img) if include_pink_right else dark
    H, W = dark.shape

    row_frac = spine_mask.mean(axis=1)
    col_frac = spine_mask.mean(axis=0)

    r_groups = _groups(np.where(row_frac > min_frac)[0], max_gap=4)
    c_groups = _groups(np.where(col_frac > min_frac)[0], max_gap=4)
    if len(r_groups) < 2 or len(c_groups) < 2:
        # Relax slightly if margins are wide
        r_groups = _groups(np.where(row_frac > 0.35)[0], max_gap=4)
        c_groups = _groups(np.where(col_frac > 0.35)[0], max_gap=4)
    if len(r_groups) < 2 or len(c_groups) < 2:
        raise ValueError(f"Could not find plot frame spines: H={r_groups}, V={c_groups}")

    top = r_groups[0]
    bot = r_groups[-1]
    left = c_groups[0]
    right = c_groups[-1]
    return {
        "r_top": top[2], "r_top_span": (top[0], top[1]),
        "r_bot": bot[2], "r_bot_span": (bot[0], bot[1]),
        "c_left": left[2], "c_left_span": (left[0], left[1]),
        "c_right": right[2], "c_right_span": (right[0], right[1]),
        "box": (left[1] + 1, top[1] + 1, right[0] - 1, bot[0] - 1),
    }


def find_left_ticks(img: np.ndarray, frame: dict, pad_lo: int = 5, pad_hi: int = 18, min_frac: float = 0.55) -> list[float]:
    """Detect outward tick-mark row centers sticking left of the left spine."""
    dark = dark_mask(img, thresh=150.0)
    c_edge = frame["c_left_span"][0]
    r0 = max(0, frame["r_top_span"][0] - 8)
    r1 = min(img.shape[0], frame["r_bot_span"][1] + 9)
    c_a = max(0, c_edge - pad_hi)
    c_b = max(c_a + 1, c_edge - pad_lo)
    strip = dark[r0:r1, c_a:c_b].mean(axis=1)
    gs = _groups(np.where(strip >= min_frac)[0] + r0, max_gap=3)
    return [g[2] for g in gs]


def find_right_ticks(img: np.ndarray, frame: dict, pad_lo: int = 5, pad_hi: int = 18,
                     min_frac: float = 0.55, use_pink: bool = True) -> list[float]:
    """Detect outward tick-mark row centers sticking right of the right spine."""
    mask = pink_mask(img) if use_pink else dark_mask(img, thresh=150.0)
    if use_pink:
        mask = mask | dark_mask(img, thresh=150.0)
    c_edge = frame["c_right_span"][1]
    r0 = max(0, frame["r_top_span"][0] - 8)
    r1 = min(img.shape[0], frame["r_bot_span"][1] + 9)
    c_a = min(img.shape[1] - 1, c_edge + pad_lo)
    c_b = min(img.shape[1], c_edge + pad_hi)
    strip = mask[r0:r1, c_a:c_b].mean(axis=1)
    gs = _groups(np.where(strip >= min_frac)[0] + r0, max_gap=3)
    return [g[2] for g in gs]


def find_bottom_ticks(img: np.ndarray, frame: dict, pad_lo: int = 5, pad_hi: int = 18, min_frac: float = 0.55) -> list[float]:
    """Detect outward tick-mark column centers sticking below the bottom spine."""
    dark = dark_mask(img, thresh=150.0)
    r_edge = frame["r_bot_span"][1]
    c0 = max(0, frame["c_left_span"][0] - 8)
    c1 = min(img.shape[1], frame["c_right_span"][1] + 9)
    r_a = min(img.shape[0] - 1, r_edge + pad_lo)
    r_b = min(img.shape[0], r_edge + pad_hi)
    strip = dark[r_a:r_b, c0:c1].mean(axis=0)
    gs = _groups(np.where(strip >= min_frac)[0] + c0, max_gap=3)
    return [g[2] for g in gs]


class LinearAxis:
    """Linear calibration between continuous pixel coordinate and printed axis value."""

    def __init__(self, marks: list[tuple[float, float]]):
        if len(marks) < 2:
            raise ValueError(f"Need >= 2 tick marks, got {marks}")
        ps = np.array([float(p) for p, _ in marks], dtype=float)
        vs = np.array([float(v) for _, v in marks], dtype=float)
        self.b, self.a = np.polyfit(ps, vs, 1)  # v = a + b * p
        pred_p = (vs - self.a) / self.b
        self.residual_px = float(np.max(np.abs(ps - pred_p)))

    def to_val(self, p: float) -> float:
        return float(self.a + self.b * float(p))

    def to_px(self, v: float) -> float:
        return float((float(v) - self.a) / self.b)

    @property
    def val_per_px(self) -> float:
        return abs(float(self.b))


def interior_mask(mask: np.ndarray, frame: dict, inset: int = 4) -> np.ndarray:
    """Restrict boolean mask to strictly inside the plot frame."""
    c0, r0, c1, r1 = frame["box"]
    out = np.zeros_like(mask, dtype=bool)
    out[r0 + inset : r1 - inset + 1, c0 + inset : c1 - inset + 1] = mask[
        r0 + inset : r1 - inset + 1, c0 + inset : c1 - inset + 1
    ]
    return out


def largest_span_component(mask: np.ndarray, min_span_frac: float = 0.5) -> np.ndarray:
    """Keep connected components with large horizontal span (isolates solid curves from labels)."""
    lab, n = ndimage.label(mask)
    if n == 0:
        return mask
    objs = ndimage.find_objects(lab)
    W = mask.shape[1]
    keep = np.zeros_like(mask, dtype=bool)
    best_span = 0
    best_idx = 0
    for idx, sl in enumerate(objs, start=1):
        if sl is None:
            continue
        span = sl[1].stop - sl[1].start
        if span > best_span:
            best_span = span
            best_idx = idx
        if span >= min_span_frac * W:
            keep[sl] |= (lab[sl] == idx)
    if not keep.any() and best_idx > 0:
        sl = objs[best_idx - 1]
        keep[sl] = (lab[sl] == best_idx)
    return keep


def _col_runs(col_bool: np.ndarray) -> list[tuple[int, int]]:
    idx = np.nonzero(col_bool)[0]
    if len(idx) == 0:
        return []
    splits = np.where(np.diff(idx) > 2)[0] + 1
    return [(int(g[0]), int(g[-1])) for g in np.split(idx, splits)]


def sample_curve_column(
    mask: np.ndarray,
    x_axis: LinearAxis,
    y_axis: LinearAxis,
    at_xs: list[float],
    half_win: int = 2,
    gap_search: int = 0,
    branch: str = "mean",
    min_half_px: float = 1.5,
    extra_px: float = 0.75,
    note_prefix: str = "",
) -> list[dict]:
    """Read y(x) at each x in `at_xs` from a cleaned curve mask.

    - `branch`: 'mean' (single curve centroid), 'top' (uppermost run), 'bot' (lowermost run),
      or a callable `(x, runs) -> (r_start, r_end)` for custom disambiguation.
    - `gap_search`: if > 0 and no ink exists in `[col - half_win, col + half_win]` (e.g. a
      dashed/dotted line gap), search up to `gap_search` px left & right and linearly interpolate
      the row center across the gap, widening the interval honestly.
    """
    H, W = mask.shape
    out = []

    def _row_at_col(c_int: int, x_val: float) -> tuple[float, float] | None:
        c0 = max(0, c_int - half_win)
        c1 = min(W, c_int + half_win + 1)
        col_any = mask[:, c0:c1].any(axis=1)
        runs = _col_runs(col_any)
        if not runs:
            return None
        if callable(branch):
            chosen = branch(x_val, runs)
            if chosen is None:
                return None
            r0, r1 = chosen
            return 0.5 * (r0 + r1) + 0.5, float(r1 - r0 + 1)
        if branch == "top":
            r0, r1 = runs[0]
            return 0.5 * (r0 + r1) + 0.5, float(r1 - r0 + 1)
        if branch == "bot":
            r0, r1 = runs[-1]
            return 0.5 * (r0 + r1) + 0.5, float(r1 - r0 + 1)
        # 'mean': pick the run closest to the weighted median, or combine if single
        if len(runs) == 1:
            r0, r1 = runs[0]
            return 0.5 * (r0 + r1) + 0.5, float(r1 - r0 + 1)
        # Multiple runs: take widest (thickest) run or median center
        r0, r1 = max(runs, key=lambda r: r[1] - r[0])
        return 0.5 * (r0 + r1) + 0.5, float(r1 - r0 + 1)

    x_err_px = 1.0 + x_axis.residual_px
    y_err_px = extra_px + y_axis.residual_px

    for x in at_xs:
        col_f = x_axis.to_px(x)
        c_int = int(round(col_f))
        if not (0 <= c_int < W):
            out.append({"x": float(x), "y": None, "found": False, "note": "x outside crop"})
            continue

        hit = _row_at_col(c_int, x)
        interp_gap = 0.0
        if hit is not None:
            row_f, thick = hit
        elif gap_search > 0:
            left_pt = None
            right_pt = None
            for d in range(half_win + 1, gap_search + 1):
                if left_pt is None and c_int - d >= 0:
                    h = _row_at_col(c_int - d, x)
                    if h is not None:
                        left_pt = (c_int - d, h[0], h[1])
                if right_pt is None and c_int + d < W:
                    h = _row_at_col(c_int + d, x)
                    if h is not None:
                        right_pt = (c_int + d, h[0], h[1])
                if left_pt is not None and right_pt is not None:
                    break
            if left_pt is not None and right_pt is not None:
                cL, rL, tL = left_pt
                cR, rR, tR = right_pt
                frac = (col_f - cL) / max(1.0, float(cR - cL))
                row_f = rL + frac * (rR - rL)
                thick = 0.5 * (tL + tR)
                interp_gap = float(cR - cL)
            elif left_pt is not None:
                row_f, thick = left_pt[1], left_pt[2]
                interp_gap = float(abs(col_f - left_pt[0])) * 2.0
            elif right_pt is not None:
                row_f, thick = right_pt[1], right_pt[2]
                interp_gap = float(abs(right_pt[0] - col_f)) * 2.0
            else:
                out.append({"x": float(x), "y": None, "found": False, "note": "no ink within gap_search"})
                continue
        else:
            out.append({"x": float(x), "y": None, "found": False, "note": "no curve ink at column"})
            continue

        # Estimate local slope |dr/dc| from +/- step columns for honest horizontal-to-vertical error
        slope = 0.0
        step = max(3, int(round(max(4.0, interp_gap * 0.5))))
        pL, cL_s = _row_at_col(max(0, c_int - step), x), max(0, c_int - step)
        pR, cR_s = _row_at_col(min(W - 1, c_int + step), x), min(W - 1, c_int + step)
        if gap_search > 0:
            if pL is None:
                for d in range(step + 1, gap_search + step + 1):
                    if c_int - d < 0:
                        break
                    h = _row_at_col(c_int - d, x)
                    if h is not None:
                        pL, cL_s = h, c_int - d
                        break
            if pR is None:
                for d in range(step + 1, gap_search + step + 1):
                    if c_int + d >= W:
                        break
                    h = _row_at_col(c_int + d, x)
                    if h is not None:
                        pR, cR_s = h, c_int + d
                        break
        if pL is not None and pR is not None and cR_s > cL_s:
            slope = abs(pR[0] - pL[0]) / float(cR_s - cL_s)
        elif pL is not None and c_int > cL_s:
            slope = abs(row_f - pL[0]) / float(c_int - cL_s)
        elif pR is not None and cR_s > c_int:
            slope = abs(pR[0] - row_f) / float(cR_s - c_int)

        half_px = max(min_half_px, 0.5 * thick + y_err_px + slope * x_err_px + 0.05 * interp_gap)
        y_val = y_axis.to_val(row_f)
        y_a = y_axis.to_val(row_f - half_px)
        y_b = y_axis.to_val(row_f + half_px)
        lo, hi = min(y_a, y_b), max(y_a, y_b)
        px_hi_y = row_f - half_px if y_axis.b < 0 else row_f + half_px
        px_lo_y = row_f + half_px if y_axis.b < 0 else row_f - half_px

        note = f"{note_prefix}thick={thick:.1f}px, slope={slope:.2f}"
        if interp_gap > 0:
            note += f", gap={interp_gap:.0f}px"
        out.append({
            "x": float(x),
            "y": y_val,
            "lo": lo,
            "hi": hi,
            "px": [float(col_f), float(row_f)],
            "px_lo": float(px_hi_y),
            "px_hi": float(px_lo_y),
            "note": note,
        })
    return out
