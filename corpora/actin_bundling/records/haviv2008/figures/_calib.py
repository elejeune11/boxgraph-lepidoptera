"""Shared frame, tick-calibration, and geometry helpers for Haviv et al. (2008) Fig. 2a-d."""
from __future__ import annotations

import math
import numpy as np
from scipy import ndimage
from boxgraph.digitize import Axis, regular_ticks

# Embedded raster panels in 2008_Haviv.pdf p. 4 are ~150 dpi (8 px per native pixel at 1200 dpi).
NATIVE_QUANT_PX = 4.0  # half a native bitmap pixel at 1200 dpi


def _group_runs(indices: np.ndarray, gap: int = 3) -> list[tuple[int, int]]:
    if len(indices) == 0:
        return []
    out: list[tuple[int, int]] = []
    s = p = int(indices[0])
    for x in indices[1:]:
        xi = int(x)
        if xi <= p + gap:
            p = xi
        else:
            out.append((s, p))
            s = p = xi
    out.append((s, p))
    return out


def find_frame(image: np.ndarray) -> dict[str, float | int]:
    """Detect the four neutral-black plot frame spines (L, R, T, B) directly from the image."""
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    gray = (R + G + B) / 3.0
    neutral = (np.abs(R - G) < 25) & (np.abs(R - B) < 25) & (np.abs(G - B) < 25)
    dark = (gray < 150) & neutral
    H, W = dark.shape

    hr = _group_runs(np.where(dark.sum(axis=1) > 0.55 * W)[0])
    vc = _group_runs(np.where(dark.sum(axis=0) > 0.55 * H)[0])
    if len(hr) < 2 or len(vc) < 2:
        raise RuntimeError(f"Could not detect 4-sided plot frame: hr={hr}, vc={vc}")
    top_r0, top_r1 = hr[0]
    bot_r0, bot_r1 = hr[-1]
    left_c0, left_c1 = vc[0]
    right_c0, right_c1 = vc[-1]

    def _refine_spine(arr_2d: np.ndarray, start: int, is_row: bool) -> tuple[float, int, int]:
        prof = np.median(arr_2d, axis=1 if is_row else 0)
        idx = np.where(prof < 150)[0]
        if len(idx) == 0:
            idx = np.where(prof < 200)[0]
        return float((idx[0] + idx[-1]) / 2.0 + 0.5 + start), int(idx[0] + start), int(idx[-1] + start)

    T, tr0, tr1 = _refine_spine(
        gray[max(0, top_r0 - 5) : top_r1 + 6, left_c1 + 20 : right_c0 - 20],
        max(0, top_r0 - 5),
        True,
    )
    Bb, br0, br1 = _refine_spine(
        gray[bot_r0 - 5 : min(H, bot_r1 + 6), left_c1 + 20 : right_c0 - 20],
        bot_r0 - 5,
        True,
    )
    L, lc0, lc1 = _refine_spine(
        gray[top_r1 + 20 : bot_r0 - 20, max(0, left_c0 - 5) : left_c1 + 6],
        max(0, left_c0 - 5),
        False,
    )
    Rr, rc0, rc1 = _refine_spine(
        gray[top_r1 + 20 : bot_r0 - 20, right_c0 - 5 : min(W, right_c1 + 6)],
        right_c0 - 5,
        False,
    )
    return {
        "L": L, "R": Rr, "T": T, "B": Bb,
        "lc0": lc0, "lc1": lc1, "rc0": rc0, "rc1": rc1,
        "tr0": tr0, "tr1": tr1, "br0": br0, "br1": br1,
    }


def scan_inward_ticks(
    image: np.ndarray,
    frame: dict[str, float | int],
    side: str,
    depth: tuple[int, int] = (4, 18),
    thresh: float = 210.0,
    max_w: int = 26,
) -> list[float]:
    """Scan inward neutral-black ticks along the bottom or left spine."""
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    gray = (R + G + B) / 3.0
    neutral = (np.abs(R - G) < 25) & (np.abs(R - B) < 25) & (np.abs(G - B) < 25)
    ngray = np.where(neutral, gray, 255.0)

    if side == "bot":
        strip = ngray[int(frame["br0"]) - depth[1] : int(frame["br0"]) - depth[0], int(frame["lc1"]) + 12 : int(frame["rc0"]) - 12]
        off = int(frame["lc1"]) + 12
        prof = strip.mean(axis=0)
    elif side == "left":
        strip = ngray[int(frame["tr1"]) + 12 : int(frame["br0"]) - 12, int(frame["lc1"]) + depth[0] : int(frame["lc1"]) + depth[1]].T
        off = int(frame["tr1"]) + 12
        prof = strip.mean(axis=0)
    else:
        raise ValueError(f"unknown side {side!r}")

    idx = np.where(prof < thresh)[0]
    if len(idx) == 0:
        return []
    out: list[float] = []
    for s, p in _group_runs(idx, gap=2):
        if p - s + 1 <= max_w:
            w = 255.0 - prof[s : p + 1]
            c = float(np.sum(np.arange(s, p + 1) * w) / np.sum(w)) + off + 0.5
            out.append(c)
    return out


def _merge_nearby(vals: list[float], tol: float = 6.0) -> list[float]:
    if not vals:
        return []
    vs = sorted(vals)
    groups = [[vs[0]]]
    for v in vs[1:]:
        if v - groups[-1][-1] <= tol:
            groups[-1].append(v)
        else:
            groups.append([v])
    return [float(np.mean(g)) for g in groups]


def calibrate_fig2a(image: np.ndarray) -> tuple[Axis, Axis, dict]:
    """Calibrate Fig. 2a: x = N_b (0, 50, ..., 400), y = P(N_b) (0.0, 0.2, ..., 1.0)."""
    f = find_frame(image)
    xt_raw = [float(f["L"])] + scan_inward_ticks(image, f, "bot", depth=(4, 18), thresh=210.0)
    xt = regular_ticks(xt_raw, tol=6.0, min_gap=180.0)
    if len(xt) < 4:
        raise RuntimeError(f"Fig. 2a expected >=4 regular x-ticks, found {xt}")
    x_axis = Axis("linear", [(p, 50.0 * i) for i, p in enumerate(xt)])

    yt_raw = scan_inward_ticks(image, f, "left", depth=(4, 18), thresh=210.0) + [float(f["B"])]
    # Adjacent 0.2 ticks have spacing ~459 px; find min gap > 300 px among detected ticks
    ys = sorted(yt_raw)
    diffs = [ys[i + 1] - ys[i] for i in range(len(ys) - 1) if ys[i + 1] - ys[i] > 300.0]
    step_y = float(min(diffs))
    ym = []
    for p in ys:
        k = round((float(f["B"]) - p) / step_y)
        if 0 <= k <= 5 and abs((float(f["B"]) - p) - k * step_y) < 12.0:
            ym.append((p, 0.2 * k))
    if len(ym) < 4:
        raise RuntimeError(f"Fig. 2a expected >=4 y-ticks, found {ym}")
    y_axis = Axis("linear", ym)
    return x_axis, y_axis, f


def calibrate_fig2b(image: np.ndarray) -> tuple[Axis, Axis, dict]:
    """Calibrate Fig. 2b: x = N_b (0..400), y = ln P(N_b) (0, -0.5, ..., -5)."""
    f = find_frame(image)
    xt_raw = [float(f["L"])] + scan_inward_ticks(image, f, "bot", depth=(4, 18), thresh=210.0)
    xt = regular_ticks(xt_raw, tol=6.0, min_gap=180.0)
    if len(xt) < 8:
        raise RuntimeError(f"Fig. 2b expected >=8 regular x-ticks, found {xt}")
    x_axis = Axis("linear", [(p, 50.0 * i) for i, p in enumerate(xt)])

    yt_raw = scan_inward_ticks(image, f, "left", depth=(4, 18), thresh=210.0)
    yt_reg = regular_ticks(yt_raw, tol=6.0, min_gap=150.0)
    if len(yt_reg) < 8:
        raise RuntimeError(f"Fig. 2b expected >=8 regular y-ticks, found {yt_reg}")
    step_y = (yt_reg[-1] - yt_reg[0]) / (len(yt_reg) - 1)
    top_y0 = min(yt_raw)
    ym = []
    for p in sorted(yt_raw):
        k = round((p - top_y0) / step_y)
        if 0 <= k <= 10 and abs((p - top_y0) - k * step_y) < 10.0:
            ym.append((p, -0.5 * k))
    if len(ym) < 9:
        raise RuntimeError(f"Fig. 2b expected >=9 y-ticks, found {ym}")
    y_axis = Axis("linear", ym)
    return x_axis, y_axis, f


def calibrate_fig2c(image: np.ndarray) -> tuple[Axis, Axis, dict]:
    """Calibrate Fig. 2c: x = phi (0..10), y = plotted rate scale (0.0..1.0)."""
    f = find_frame(image)
    xt_raw = [float(f["L"])] + scan_inward_ticks(image, f, "bot", depth=(4, 18), thresh=215.0) + [float(f["R"])]
    xt = regular_ticks(xt_raw, tol=6.0, min_gap=180.0)
    if len(xt) != 11:
        raise RuntimeError(f"Fig. 2c expected 11 x-ticks (0..10), found {len(xt)}: {xt}")
    x_axis = Axis("linear", [(p, float(i)) for i, p in enumerate(xt)])

    yt_a = scan_inward_ticks(image, f, "left", depth=(3, 12), thresh=215.0, max_w=24)
    yt_b = scan_inward_ticks(image, f, "left", depth=(6, 20), thresh=215.0, max_w=24)
    yt_all = _merge_nearby(yt_a + yt_b + [float(f["B"])], tol=6.0)
    # Regular step for 0.2 units is ~455 px (between 1.0, 0.8, 0.6, 0.4, 0.2, 0.0)
    # Find the pair spacing closest to (B - T) / 5.5
    approx_step = (float(f["B"]) - float(f["T"])) / 5.54
    ym = []
    for p in yt_all:
        k = round((float(f["B"]) - p) / approx_step)
        if 0 <= k <= 5 and abs((float(f["B"]) - p) - k * approx_step) < 20.0:
            ym.append((p, 0.2 * k))
    # Refine step from matched ticks and drop any outlier > 6 px
    if len(ym) < 4:
        raise RuntimeError(f"Fig. 2c expected >=4 y-ticks, found {ym} from {yt_all}")
    y_axis = Axis("linear", ym)
    if y_axis.residual_px > 5.0:
        ym = [(p, v) for p, v in ym if abs(p - y_axis.to_pixel(v)) <= 5.0]
        y_axis = Axis("linear", ym)
    return x_axis, y_axis, f


def calibrate_fig2d(image: np.ndarray) -> tuple[Axis, Axis, dict]:
    """Calibrate Fig. 2d: x = phi (0..10), y = N_peak (0, 5, ..., 50)."""
    f = find_frame(image)
    xt_raw = [float(f["L"])] + scan_inward_ticks(image, f, "bot", depth=(4, 18), thresh=215.0) + [float(f["R"])]
    xt = regular_ticks(xt_raw, tol=6.0, min_gap=180.0)
    if len(xt) != 11:
        raise RuntimeError(f"Fig. 2d expected 11 x-ticks (0..10), found {len(xt)}: {xt}")
    x_axis = Axis("linear", [(p, float(i)) for i, p in enumerate(xt)])

    yt_raw = scan_inward_ticks(image, f, "left", depth=(4, 16), thresh=215.0, max_w=24)
    yt_reg = regular_ticks(yt_raw, tol=6.0, min_gap=180.0)
    if len(yt_reg) < 8:
        raise RuntimeError(f"Fig. 2d expected >=8 regular y-ticks, found {yt_reg}")
    # Top tick is 50, descending by 5 per step down to 10 (plus bottom spine B = 0)
    step_y = (yt_reg[-1] - yt_reg[0]) / (len(yt_reg) - 1)
    ym = [(p, 50.0 - 5.0 * round((p - yt_reg[0]) / step_y)) for p in yt_reg]
    if abs((float(f["B"]) - yt_reg[0]) - 10.0 * step_y) < 15.0:
        ym.append((float(f["B"]), 0.0))
    y_axis = Axis("linear", ym)
    return x_axis, y_axis, f


def value_interval(axis: Axis, px: float, half_px: float) -> tuple[float, float, float]:
    """Convert pixel coordinate +/- half_px into (val, lo, hi) with lo <= val <= hi."""
    val = float(axis.to_value(px))
    a = float(axis.to_value(px - half_px))
    b = float(axis.to_value(px + half_px))
    return val, min(a, b), max(a, b)


def col_runs(mask_col: np.ndarray, row_offset: int) -> list[tuple[float, float, int, int]]:
    """Return list of (center_row, half_height, r0, r1) for True runs in a 1D boolean column."""
    idx = np.where(mask_col)[0]
    if len(idx) == 0:
        return []
    out = []
    for s, p in _group_runs(idx, gap=2):
        r0, r1 = s + row_offset, p + row_offset
        out.append(((r0 + r1) / 2.0 + 0.5, (r1 - r0 + 1) / 2.0, r0, r1))
    return out


def detect_fig2c_marker_columns(image: np.ndarray, f: dict) -> list[float]:
    """Detect the 5 marker columns in Fig. 2c from the 5 black hollow circles (fast 1/4-res)."""
    sub = image[::4, ::4].astype(float)
    sR, sG, sB = sub[:, :, 0], sub[:, :, 1], sub[:, :, 2]
    sgray = (sR + sG + sB) / 3.0
    sint = np.zeros(sgray.shape, dtype=bool)
    sint[int(f["tr1"]) // 4 + 3 : int(f["br0"]) // 4 - 3, int(f["lc1"]) // 4 + 3 : int(f["rc0"]) // 4 - 3] = True
    bk = (sgray < 160) & (sR - sG < 25) & sint
    bk_filled = ndimage.binary_fill_holes(bk)
    yy, xx = np.ogrid[-4:5, -4:5]
    disk9 = (xx * xx + yy * yy) <= 17
    bk_circ = ndimage.binary_opening(bk_filled, structure=disk9)
    lbl, _ = ndimage.label(bk_circ)
    cols = []
    for i, sl in enumerate(ndimage.find_objects(lbl), 1):
        m = lbl[sl] == i
        _, cx = ndimage.center_of_mass(m)
        cols.append(float((cx + sl[1].start) * 4.0 + 2.0))
    cols.sort()
    return cols


def detect_fig2d_marker_columns(image: np.ndarray, f: dict) -> list[float]:
    """Detect the 5 marker columns in Fig. 2d from the 5 grey-filled black squares (fast 1/4-res)."""
    sub = image[::4, ::4].astype(float)
    sgray = sub.mean(axis=2)
    sint = np.zeros(sgray.shape, dtype=bool)
    sint[int(f["tr1"]) // 4 + 3 : int(f["br0"]) // 4 - 3, int(f["lc1"]) // 4 + 3 : int(f["rc0"]) // 4 - 3] = True
    sq8 = np.ones((8, 8), dtype=bool)
    d_sq = ndimage.binary_opening((sgray < 185) & sint, structure=sq8)
    lbl, _ = ndimage.label(d_sq)
    cols = []
    for i, sl in enumerate(ndimage.find_objects(lbl), 1):
        m = lbl[sl] == i
        _, cx = ndimage.center_of_mass(m)
        cols.append(float((cx + sl[1].start) * 4.0 + 2.0))
    cols.sort()
    return cols


def find_hline_row_red(image: np.ndarray, f: dict) -> tuple[float, int, int]:
    """Detect the solid red horizontal line in Fig. 2c directly from row projection."""
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    red = (R - G > 20) & (R - B > 20)
    frac = red[:, int(f["lc1"]) + 20 : int(f["rc0"]) - 20].mean(axis=1)
    rows = np.where(frac > 0.70)[0]
    if len(rows) == 0:
        raise RuntimeError("Could not detect solid red horizontal line in Fig. 2c")
    return float((rows[0] + rows[-1]) / 2.0 + 0.5), int(rows[0]), int(rows[-1])
