"""Calibrated image reader for S10 Figure (p. 35 / SI-16): Rg(t)/Rg(0) time series at alpha:A = 0.1
across 6 M:A ratios (0.0225, 0.045, 0.09, 0.18, 0.225, 0.675) for (a) Unipolar and (b) Apolar
initial configurations (900 dpi).

Uses orthogonal RGB channel separation (black, red-including-green-tinted-rust, bright green)
accounting for alpha-blended translucent SD bands, local linear window fitting to bridge dashed
segments with sub-pixel accuracy, and direct SD band edge measurement where bands are unobstructed.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _common as C


def calibrate_s10(image: np.ndarray, polarity: str):
    info = C.find_matlab_box(
        image, has_top=False, has_right=False, dark_thresh=150.0,
        min_tick_len=5, max_tick_len=55,
    )
    sp = info["spines"]
    xt = [sp["left"]] + info["x_ticks_bot"]
    xa = C.fit_axis(xt, [0.0, 500.0, 1000.0, 1500.0, 2000.0], max_residual_px=3.0)
    yvals_all = [1.4, 1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7] if polarity == "unipolar" else [2.0, 1.8, 1.6, 1.4, 1.2, 1.0, 0.8]
    y_top = info["y_ticks_left"][0]
    y_bot = sp["bottom"]
    expected = np.linspace(y_top, y_bot, len(yvals_all))
    all_cands = info["y_ticks_left"] + [y_bot]
    mpx, mv = [], []
    for ey, val in zip(expected, yvals_all):
        near = [t for t in all_cands if abs(t - ey) <= 10.0]
        if near:
            mpx.append(min(near, key=lambda t: abs(t - ey)))
            mv.append(val)
    ya = C.fit_axis(mpx, mv, max_residual_px=3.0)
    return info, xa, ya


def _apply_interior_mask(m: np.ndarray, box: tuple[int, int, int, int], mask_title: bool = False) -> np.ndarray:
    c0, r0, c1, r1 = box
    W = c1 - c0
    H = r1 - r0
    out = np.zeros_like(m)
    bot_pad = 44 if mask_title else 8
    left_pad = 34 if mask_title else 8
    out[r0 + 8:r1 - bot_pad, c0 + left_pad:c1 + 8] = m[r0 + 8:r1 - bot_pad, c0 + left_pad:c1 + 8]
    # Mask legend box (upper left)
    out[r0:r0 + int(0.42 * H), c0:c0 + int(0.28 * W)] = False
    if mask_title:
        # Panel title text is pure black ("a) Unipolar..." / "b) Apolar...", upper center)
        out[r0:r0 + int(0.24 * H), c0 + int(0.25 * W):c0 + int(0.76 * W)] = False
    return out


def _color_line_mask(image: np.ndarray, color: str, box: tuple[int, int, int, int]) -> np.ndarray:
    r = image[:, :, 0].astype(np.int16)
    g = image[:, :, 1].astype(np.int16)
    b = image[:, :, 2].astype(np.int16)
    if color == "black":
        m = (r < 95) & (g < 115) & (b < 75) & ((r - g) < 40)
    elif color == "red":
        m = (r > 120) & ((r - g) > 65) & (b < 80)
    elif color == "green":
        m = (g > 180) & (b < 85) & ((g - b) > 130) & ((g - r) > 65)
    else:
        raise ValueError(f"Unknown color {color}")
    return _apply_interior_mask(m, box, mask_title=(color == "black"))


def _band_mask(image: np.ndarray, color: str, box: tuple[int, int, int, int]) -> np.ndarray:
    """Shaded translucent SD band mask (including the line itself and mild alpha-blend overlap)."""
    r = image[:, :, 0].astype(np.int16)
    g = image[:, :, 1].astype(np.int16)
    b = image[:, :, 2].astype(np.int16)
    if color == "green":
        # Light green SD band (even where slightly overlaid with pink SD band)
        m = (g > 165) & ((g - b) > 35) & (b < 195)
    elif color == "red":
        m = (r > 120) & ((r - g) > 25) & ((r - b) > 5)
    elif color == "black":
        lum = r + g + b
        sat = np.maximum(np.maximum(r, g), b) - np.minimum(np.minimum(r, g), b)
        m = (lum < 695) & (sat < 26)
    else:
        raise ValueError(color)
    return _apply_interior_mask(m, box, mask_title=(color == "black"))


def _cluster_runs_in_window(mask: np.ndarray, col_c: int, half_w: int,
                            r_lo: int, r_hi: int, gap_r: int = 10) -> list[tuple[np.ndarray, np.ndarray]]:
    H, W = mask.shape
    c_a = max(0, col_c - half_w)
    c_b = min(W, col_c + half_w + 1)
    pts_c, pts_r = [], []
    for c in range(c_a, c_b):
        rs = np.flatnonzero(mask[r_lo:r_hi, c]) + r_lo
        if rs.size == 0:
            continue
        sp = np.flatnonzero(np.diff(rs) > 3)
        st = np.concatenate(([0], sp + 1))
        en = np.concatenate((sp, [rs.size - 1]))
        for a, b in zip(st, en):
            h = rs[b] - rs[a] + 1
            if h <= 30:
                pts_c.append(c + 0.5)
                pts_r.append(0.5 * (rs[a] + rs[b]) + 0.5)
    if not pts_c:
        return []
    order = np.argsort(pts_c)
    tracks: list[list[tuple[float, float]]] = []
    for idx in order:
        cx, ry = pts_c[idx], pts_r[idx]
        best_t, best_d = None, gap_r + 1.0
        for ti, tr in enumerate(tracks):
            dx = cx - tr[-1][0]
            if dx > 24.0:
                continue
            if len(tr) >= 2:
                xs = np.array([p[0] for p in tr[-5:]])
                ys = np.array([p[1] for p in tr[-5:]])
                if xs[-1] > xs[0]:
                    m = (ys[-1] - ys[0]) / (xs[-1] - xs[0])
                    pred = ys[-1] + m * (cx - xs[-1])
                else:
                    pred = ys[-1]
            else:
                pred = tr[-1][1]
            d = abs(ry - pred)
            if d < best_d:
                best_d = d
                best_t = ti
        if best_t is not None:
            tracks[best_t].append((cx, ry))
        else:
            tracks.append([(cx, ry)])
    clusters = []
    for tr in tracks:
        if len(tr) >= 4:
            cs = np.array([p[0] for p in tr], dtype=float)
            rs_arr = np.array([p[1] for p in tr], dtype=float)
            clusters.append((cs, rs_arr))

    def eval_at(cl):
        cs, rs_arr = cl
        if cs.max() - cs.min() >= 3.0:
            m, b = np.polyfit(cs, rs_arr, 1)
            return float(m * col_c + b)
        return float(np.median(rs_arr))

    clusters.sort(key=eval_at)
    return clusters


def _measure_sd_band(band: np.ndarray, col_c: int, row_c: float, r_lo: int, r_hi: int,
                     max_half: int = 320, mirror_side: str | None = None) -> tuple[float | None, float | None]:
    H, W = band.shape
    c_a = max(0, col_c - 2)
    c_b = min(W, col_c + 3)
    col_has = band[:, c_a:c_b].mean(axis=1) >= 0.4
    r0_i = int(round(row_c))
    rt = r0_i
    gap = 0
    for r in range(r0_i, max(r_lo, r0_i - max_half) - 1, -1):
        if col_has[r]:
            rt = r
            gap = 0
        else:
            gap += 1
            if gap > 8:
                break
    rb = r0_i
    gap = 0
    for r in range(r0_i, min(r_hi, r0_i + max_half)):
        if col_has[r]:
            rb = r
            gap = 0
        else:
            gap += 1
            if gap > 8:
                break
    if mirror_side == "upper":
        d = row_c - rt
        return float(rt), float(row_c + d)
    if mirror_side == "lower":
        d = rb - row_c
        return float(row_c - d), float(rb)
    return float(rt), float(rb)


def read_s10_series(image: np.ndarray, polarity: str, ma: float, at: list[float],
                    read_sd: bool = False) -> list[dict]:
    """Read one M:A time series from S10(a) unipolar or S10(b) apolar."""
    info, xa, ya = calibrate_s10(image, polarity)
    c0, r0, c1, r1 = info["box"]

    color_map = {
        0.0225: ("black", "solid"),
        0.045: ("black", "dashed"),
        0.09: ("red", "solid"),
        0.18: ("red", "dashed"),
        0.225: ("green", "solid"),
        0.675: ("green", "dashed"),
    }
    color, style = color_map[ma]
    mask = _color_line_mask(image, color, info["box"])
    band = _band_mask(image, color, info["box"]) if read_sd else None

    clamp = (0.7, 1.4) if polarity == "unipolar" else (0.8, 2.0)
    half_w = 24 if style == "dashed" else 14

    bot_lim = r1 - (44 if color == "black" else 8)
    out = []
    for t in at:
        pc = xa.to_pixel(t)
        col_c = int(round(pc))
        hw = half_w
        if col_c + hw > c1 + 4:
            c_win = col_c - hw
            clusters = _cluster_runs_in_window(mask, c_win, hw, r0 + 8, bot_lim, gap_r=10)
        else:
            clusters = _cluster_runs_in_window(mask, col_c, hw, r0 + 8, bot_lim, gap_r=10)

        if not clusters:
            out.append({"x": float(t), "y": None, "note": f"S10 {polarity} M:A={ma}: no {color} curve at t={t}"})
            continue

        if len(clusters) == 1:
            chosen = clusters[0]
        else:
            if polarity == "unipolar":
                if ma in (0.0225, 0.09):
                    chosen = clusters[0]
                elif ma in (0.045, 0.18):
                    chosen = clusters[-1]
                elif ma == 0.225:
                    chosen = clusters[0] if t < 1000.0 else clusters[-1]
                else:
                    chosen = clusters[-1] if t < 1000.0 else clusters[0]
            else:
                if ma in (0.045, 0.18, 0.225):
                    chosen = clusters[0]
                else:
                    chosen = clusters[-1]

        cs, rs_arr = chosen
        if cs.max() - cs.min() >= 4.0:
            m, b = np.polyfit(cs, rs_arr, 1)
            row_val = float(m * pc + b)
            slope = float(m)
            res = float(np.median(np.abs(rs_arr - (m * cs + b))))
        else:
            row_val = float(np.median(rs_arr))
            slope = 0.0
            res = 2.0

        rt, rb = None, None
        sd_note = ""
        if read_sd and band is not None:
            mirror = None
            skip_sd = False
            if polarity == "apolar" and ma == 0.09:
                if 1450.0 <= t <= 1650.0:
                    mirror = "upper"
                elif t > 1650.0:
                    skip_sd = True
            elif polarity == "apolar" and ma == 0.675 and t >= 1500.0:
                mirror = "upper"
            if not skip_sd:
                rt, rb = _measure_sd_band(band, min(col_c, c1 - 4), row_val, r0 + 8, bot_lim, mirror_side=mirror)
                if rt is not None and rb is not None:
                    sd_note = "; SD band" + (f" ({mirror} mirrored)" if mirror else "")

        note = f"S10 {polarity} alpha:A=0.1 M:A={ma:g} ({color} {style}{sd_note}; rx={xa.residual_px:.2f}px, ry={ya.residual_px:.2f}px)"
        out.append(
            C.make_reading(
                t, row_val, ya, pc, half_px=max(2.5, res + 1.5), x_axis=xa,
                slope_px=slope, note=note, bar_top_px=rt, bar_bot_px=rb,
                bar_err_px=3.0, clamp_y=clamp,
            )
        )
    return out
