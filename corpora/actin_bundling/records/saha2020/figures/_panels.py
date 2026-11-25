"""Panel-specific calibration and extraction routines for saha2020 figures."""

from __future__ import annotations

import math
import numpy as np

from _common import (
    LinearAxis,
    _group_indices,
    column_runs,
    estimate_slope,
    find_frame,
    inward_ticks,
    make_reading,
    read_branch_from_mask,
    slice_ticks,
)


# ---------------------------------------------------------------------------
# Figure 1 (d, g, l): 1D concentration line scans phi_1(r), phi_2(r)
# ---------------------------------------------------------------------------

def read_fig1(image: np.ndarray, at: list[float], panel: str, species: str) -> list[dict]:
    arr = image.astype(np.float32)
    r, g, b = arr[:, :, 0], arr[:, :, 1], arr[:, :, 2]
    gray = arr.mean(axis=2)
    sat = arr.max(axis=2) - arr.min(axis=2)

    c_l, r_t, c_r, r_b = find_frame(image, thr=205.0, min_frac=0.45, neutral_only=False)
    ink = (gray < 205.0) & (sat < 30.0)

    c0, r0, c1, r1 = c_l[2], r_t[2], c_r[2], r_b[2]
    bot_t = slice_ticks(ink, r_b[0] - 15, c_l[1] + 15, c_r[0] - 15, is_x=True, max_width=28)
    left_t = slice_ticks(ink, c_l[1] + 15, r_t[1] + 15, r_b[0] - 15, is_x=False, max_width=28)
    right_t = slice_ticks(ink, c_r[0] - 15, r_t[1] + 15, r_b[0] - 15, is_x=False, max_width=28)

    if panel == "d":
        if len(bot_t) < 3:
            raise RuntimeError(f"fig1d: expected 3 interior x ticks, got {bot_t}")
        xa = LinearAxis([(c0, -1.0), (bot_t[0], -0.5), (bot_t[1], 0.0), (bot_t[2], 0.5), (c1, 1.0)])
        # Filter y ticks to the 3 regularly spaced ticks (0.6, 0.4, 0.2) plus bottom spine (0.0)
        # On right spine, ticks are cleanly separated from blue curve at r=+1 (row ~253 vs 0.6 at ~299)
        y_cands = sorted(t for t in left_t if t > r_t[1] + 120)
        if len(y_cands) < 3:
            y_cands = sorted(t for t in right_t if t > r_t[1] + 120)
        ya = LinearAxis([(y_cands[0], 0.6), (y_cands[1], 0.4), (y_cands[2], 0.2), (r1, 0.0)])
    elif panel == "g":
        if len(bot_t) < 4:
            raise RuntimeError(f"fig1g: expected 4 x ticks, got {bot_t}")
        xa = LinearAxis([(bot_t[0], -0.8), (bot_t[1], -0.6), (bot_t[2], -0.4), (bot_t[3], -0.2)])
        # Right spine has crisp ticks at 0.6, 0.4, 0.2 (~314 px spacing from bottom spine r1 = 0.0)
        y_all = sorted(set([round(t, 1) for t in right_t + left_t]))
        best = None
        for y02 in y_all:
            st = r1 - y02
            if not (0.24 * (r1 - r0) < st < 0.34 * (r1 - r0)):
                continue
            y04 = min(y_all, key=lambda t: abs(t - (r1 - 2.0 * st)))
            y06 = min(y_all, key=lambda t: abs(t - (r1 - 3.0 * st)))
            err = abs(y04 - (r1 - 2.0 * st)) + abs(y06 - (r1 - 3.0 * st))
            if best is None or err < best[0]:
                best = (err, (y06, y04, y02))
        y3 = best[1]
        ya = LinearAxis([(y3[0], 0.6), (y3[1], 0.4), (y3[2], 0.2)])
    elif panel == "l":
        if len(bot_t) < 3:
            raise RuntimeError(f"fig1l: expected 3 interior x ticks, got {bot_t}")
        xa = LinearAxis([(c0, -1.0), (bot_t[0], -0.5), (bot_t[1], 0.0), (bot_t[2], 0.5), (c1, 1.0)])
        # 3 y ticks (0.5, 0.3, 0.1) appear on both left and right spines at identical rows
        matched = []
        for lt in left_t:
            if any(abs(lt - rt) < 5.0 for rt in right_t):
                matched.append(lt)
        # Keep the 3 ticks spaced by ~445 px (matching ~0.40 of box height)
        H_box = r1 - r0
        best = None
        for i in range(len(matched)):
            for j in range(i + 1, len(matched)):
                st = matched[j] - matched[i]
                if not (0.30 * H_box < st < 0.50 * H_box):
                    continue
                for k in range(j + 1, len(matched)):
                    err = abs((matched[k] - matched[j]) - st)
                    if best is None or err < best[0]:
                        best = (err, (matched[i], matched[j], matched[k]))
        y3 = best[1]
        ya = LinearAxis([(y3[0], 0.5), (y3[1], 0.3), (y3[2], 0.1)])
    else:
        raise ValueError(f"Unknown fig1 panel {panel}")

    if species == "phi1":
        # Yellow-green curve
        mask = (g - b > 28.0) & (g > 130.0) & (r > 90.0)
    else:
        # Blue curve
        mask = (b - r > 35.0) & (b > 130.0) & (r < 175.0)

    # Mask out interior legend box so legend swatches never interfere
    H_box, W_box = r1 - r0, c1 - c0
    if panel in ("d", "g"):
        mask[int(r0):int(r0 + 0.18 * H_box), int(c0 + 0.20 * W_box):int(c0 + 0.62 * W_box)] = False
    elif panel == "l":
        mask[int(r1 - 0.18 * H_box):int(r1), int(c0 + 0.22 * W_box):int(c0 + 0.56 * W_box)] = False

    return read_branch_from_mask(
        mask, xa, ya, at, int(r_t[1] + 2), int(r_b[0] - 2),
        pick="thickest", search=4, note_prefix=f"Fig1({panel}) {species}"
    )


# ---------------------------------------------------------------------------
# Figure 2 (a, b, c): Spinodal C_r and complex-eigenvalue C_i boundaries
# ---------------------------------------------------------------------------

def _calibrate_fig2(image: np.ndarray):
    arr = image.astype(np.float32)
    gray = arr.mean(axis=2)
    sat = arr.max(axis=2) - arr.min(axis=2)
    c_l, r_t, c_r, r_b = find_frame(image, thr=160.0, min_frac=0.45)
    ink = (gray < 160.0) & (sat < 40.0)
    bt = inward_ticks(ink, r_b[0], -1, c_l[1] + 15, c_r[0] - 15, is_x=True, min_len=20, max_len=55)
    lt = inward_ticks(ink, c_l[1], +1, r_t[1] + 15, r_b[0] - 15, is_x=False, min_len=20, max_len=55)
    if len(bt) != 3 or len(lt) != 3:
        raise RuntimeError(f"fig2 calibration failed: bt={bt}, lt={lt}")
    xa = LinearAxis([(bt[0][0], -0.5), (bt[1][0], 0.0), (bt[2][0], 0.5)])
    ya = LinearAxis([(lt[0][0], 0.5), (lt[1][0], 0.0), (lt[2][0], -0.5)])
    return xa, ya, (c_l, r_t, c_r, r_b)


def _read_shaded_boundary(fill_mask: np.ndarray, xa: LinearAxis, ya: LinearAxis,
                          at: list[float], r_lo: int, r_hi: int, side: str,
                          max_gap: int = 20, min_span: int = 30,
                          stroke_half_px: float = 2.0, note_prefix: str = "") -> list[dict]:
    """Read upper or lower boundary of a shaded spinodal interior (bridging internal tie lines)."""
    out = []
    for xv in at:
        pc = xa.to_pixel(xv)
        col = int(round(pc))
        # Find longest contiguous run of shaded interior in [r_lo, r_hi]
        idx = np.where(fill_mask[r_lo:r_hi, col])[0] + r_lo
        gs = _group_indices(idx, max_gap=max_gap)
        valid = [g for g in gs if (g[1] - g[0]) >= min_span]
        if not valid:
            out.append({"x": float(xv), "y": None, "note": f"no shaded interior near x={xv:g}"})
            continue
        best = max(valid, key=lambda g: g[1] - g[0])
        edge_row = float(best[0]) - 1.0 if side == "upper" else float(best[1]) + 1.0

        # Estimate local slope from col +/- 4
        def _edge_at(c_off):
            ix = np.where(fill_mask[r_lo:r_hi, col + c_off])[0] + r_lo
            gg = [g for g in _group_indices(ix, max_gap=max_gap) if (g[1] - g[0]) >= min_span]
            if not gg:
                return None
            b = max(gg, key=lambda g: g[1] - g[0])
            return float(b[0]) - 1.0 if side == "upper" else float(b[1]) + 1.0

        el, er = _edge_at(-4), _edge_at(4)
        slope = (er - el) / 8.0 if (el is not None and er is not None) else 0.0
        note = f"{note_prefix}; cal_misfit=({xa.residual_px:.2f},{ya.residual_px:.2f})px"
        out.append(make_reading(xv, pc, edge_row, stroke_half_px, slope, xa, ya, extra_px=1.5, note=note))
    return out


def read_fig2(image: np.ndarray, at: list[float], panel: str, branch: str) -> list[dict]:
    xa, ya, (c_l, r_t, c_r, r_b) = _calibrate_fig2(image)
    arr = image.astype(np.float32)
    r, g, b = arr[:, :, 0], arr[:, :, 1], arr[:, :, 2]
    r_lo = int(r_t[1] + 5)
    r_hi = int(r_b[0] - 5)

    if panel == "a":
        if branch in ("cr_upper", "cr_lower"):
            # Pale blue shaded spinodal interior (excluding white binodal region)
            fill = (b - r >= 6.0) & (r >= 210.0) & (r <= 246.0) & (b >= 230.0)
            side = "upper" if branch == "cr_upper" else "lower"
            return _read_shaded_boundary(fill, xa, ya, at, r_lo, r_hi, side, note_prefix=f"Fig2(a) {branch}")
        elif branch == "ci_lower":
            orange = (r > 170.0) & (r - b > 35.0) & (g < 200.0)
            r_mid = int(round(ya.to_pixel(0.0)))
            return read_branch_from_mask(
                orange, xa, ya, at, r_lo, r_mid, pick="bottom", search=3,
                extra_px=1.0, note_prefix="Fig2(a) C_i lower"
            )
    elif panel == "b":
        if branch == "cr_ul_upper":
            fill = (g - b >= 6.0) & (g >= r + 1.0) & (r >= 215.0) & (g >= 230.0) & (b <= 242.0)
            r_mid = int(round(ya.to_pixel(0.0)))
            return _read_shaded_boundary(fill, xa, ya, at, r_lo, r_mid, "upper", note_prefix="Fig2(b) C_r UL upper")
        elif branch == "cr_lr_lower":
            fill = (g - b >= 6.0) & (g >= r + 1.0) & (r >= 215.0) & (g >= 230.0) & (b <= 242.0)
            r_top_lr = int(round(ya.to_pixel(0.15)))
            return _read_shaded_boundary(fill, xa, ya, at, r_top_lr, r_hi, "lower", note_prefix="Fig2(b) C_r LR lower")
        elif branch == "ci_ll_upper":
            orange = (r > 160.0) & (r - b > 50.0) & (g < 170.0)
            r_mid = int(round(ya.to_pixel(0.0)))
            return read_branch_from_mask(
                orange, xa, ya, at, r_mid, r_hi, pick="top", search=3,
                extra_px=1.0, note_prefix="Fig2(b) C_i LL upper"
            )
    elif panel == "c":
        if branch == "cr_alpha0_upper":
            # Outermost purple contour (alpha=0): highest non-white colored curve in column
            colored = (arr.max(axis=2) - arr.min(axis=2) > 25.0) & (arr.min(axis=2) < 220.0)
            # Exclude legend on bottom-left (phi1 < -0.35 or phi2 < 0.08 when phi1 < -0.1)
            r_leg = int(round(ya.to_pixel(0.10)))
            c_leg = int(round(xa.to_pixel(-0.15)))
            colored[r_leg:, :c_leg] = False
            return read_branch_from_mask(
                colored, xa, ya, at, r_lo, r_hi, pick="top", search=3,
                extra_px=1.0, note_prefix="Fig2(c) alpha=0 C_r upper"
            )
        elif branch in ("cr_alpha09_upper", "cr_alpha09_lower"):
            red = (r > 190.0) & (g < 90.0) & (b < 80.0)
            r_mid = int(round(ya.to_pixel(0.15)))
            pick = "top" if branch == "cr_alpha09_upper" else "bottom"
            return read_branch_from_mask(
                red, xa, ya, at, r_lo, r_mid, pick=pick, search=3,
                extra_px=1.0, note_prefix=f"Fig2(c) alpha=0.9 {pick}"
            )
    raise ValueError(f"Unknown fig2 panel={panel} branch={branch}")


# ---------------------------------------------------------------------------
# Figure 3 (a, c): Coarsening R(t) and polar-order flux J_i(t)
# ---------------------------------------------------------------------------

def read_fig3a(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    arr = image.astype(np.float32)
    r, g, b = arr[:, :, 0], arr[:, :, 1], arr[:, :, 2]
    gray = arr.mean(axis=2)
    sat = arr.max(axis=2) - arr.min(axis=2)
    c_l, r_t, c_r, r_b = find_frame(image, thr=160.0, min_frac=0.45)
    ink = (gray < 160.0) & (sat < 35.0)

    bt = inward_ticks(ink, r_b[0], -1, c_l[1] + 15, c_r[0] - 15, is_x=True, min_len=35, max_len=80)
    lt = inward_ticks(ink, c_l[1], +1, r_t[1] + 15, r_b[0] - 15, is_x=False, min_len=15, max_len=80)
    if len(bt) < 3:
        raise RuntimeError(f"fig3a x major ticks failed: {bt}")
    # Major x ticks: 10^-1, 10^0, 10^1, plus right spine 10^2
    xa = LinearAxis([(bt[0][0], 0.1), (bt[1][0], 1.0), (bt[2][0], 10.0), (c_r[2], 100.0)], scale="log10")
    # Y ticks: longest interior tick is 1.0; above it is 2.0; below are 0.9..0.2; bottom spine is 0.1
    y_major = max(lt, key=lambda t: t[1])[0]
    y_below = sorted(t[0] for t in lt if t[0] > y_major + 15)
    y_above = sorted(t[0] for t in lt if t[0] < y_major - 15)
    y_marks = [(y_major, 1.0), (r_b[2], 0.1)]
    if y_above:
        y_marks.append((y_above[-1], 2.0))
    if len(y_below) == 8:
        for pos, val in zip(y_below, [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2]):
            y_marks.append((pos, val))
    ya = LinearAxis(y_marks, scale="log10")

    if series == "alpha0":
        mask = (b - r > 55.0) & (b > 110.0) & (r < 120.0)
    else:
        mask = (r - b > 65.0) & (r > 130.0) & (g < 100.0)
    # Exclude legend in top-left and red letters A-E near bottom
    H_box, W_box = r_b[0] - r_t[1], c_r[0] - c_l[1]
    mask[int(r_t[1]):int(r_t[1] + 0.30 * H_box), int(c_l[1]):int(c_l[1] + 0.35 * W_box)] = False
    mask[int(r_b[0] - 0.16 * H_box):int(r_b[0]), :] = False

    return read_branch_from_mask(
        mask, xa, ya, at, int(r_t[1] + 5), int(r_b[0] - 5),
        pick="thickest", search=6, max_gap=4, extra_px=1.0, note_prefix=f"Fig3(a) {series}"
    )


def read_fig3c(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    arr = image.astype(np.float32)
    r, g, b = arr[:, :, 0], arr[:, :, 1], arr[:, :, 2]
    gray = arr.mean(axis=2)
    sat = arr.max(axis=2) - arr.min(axis=2)
    c_l, r_t, c_r, r_b = find_frame(image, thr=180.0, min_frac=0.45, neutral_only=False)
    ink = (gray < 190.0) & (sat < 45.0)

    # Top spine inward major ticks (longer than minor ticks -> sample at +40px below top spine!)
    major_x = slice_ticks(ink, r_t[1] + 40, c_l[1] + 15, c_r[0] - 15, is_x=True, max_width=28)
    if len(major_x) != 3:
        raise RuntimeError(f"fig3c expected 3 major x ticks (0.1, 1, 10), got {major_x}")
    xa = LinearAxis([(major_x[0], 0.1), (major_x[1], 1.0), (major_x[2], 10.0)], scale="log10")
    # Y axis: 0.2 at top spine, 0.0 at bottom spine (where blue alpha=0.15 curve sits)
    ya = LinearAxis([(r_t[2], 0.2), (r_b[2], 0.0)], scale="linear")

    if series == "Jx_095":
        # Dark red curve (R ~140-205, G < 50, B < 50)
        mask = (r > 110.0) & (r < 215.0) & (r - g > 70.0) & (g < 65.0) & (b < 65.0)
        # Fallback to any red if completely overlapping bright red
        red_any = (r > 130.0) & (r - g > 75.0) & (g < 70.0) & (b < 70.0)
    elif series == "Jy_095":
        # Bright red curve (R >= 210, G < 60, B < 60)
        mask = (r >= 205.0) & (r - g > 120.0) & (g < 75.0) & (b < 75.0)
        red_any = (r > 130.0) & (r - g > 75.0) & (g < 70.0) & (b < 70.0)
    elif series == "J_015":
        # Blue curves sitting along bottom spine (y = 0)
        mask = (b > 130.0) & (b - r > 60.0)
        red_any = mask
    else:
        raise ValueError(f"Unknown fig3c series {series}")

    # Exclude legend box in upper left
    H_box, W_box = r_b[0] - r_t[1], c_r[0] - c_l[1]
    mask[int(r_t[1]):int(r_t[1] + 0.55 * H_box), int(c_l[1]):int(c_l[1] + 0.45 * W_box)] = False
    red_any[int(r_t[1]):int(r_t[1] + 0.55 * H_box), int(c_l[1]):int(c_l[1] + 0.45 * W_box)] = False

    r_lo = int(r_t[1] + 3)
    r_hi = int(r_b[1] + 3) if series == "J_015" else int(r_b[0] + 5)

    out = []
    for xv in at:
        pc = xa.to_pixel(xv)
        col = int(round(pc))
        runs = column_runs(mask, col, r_lo, r_hi, max_gap=4)
        used_mask = mask
        note_extra = ""
        if not runs and series in ("Jx_095", "Jy_095"):
            runs = column_runs(red_any, col, r_lo, r_hi, max_gap=4)
            used_mask = red_any
            note_extra = "overlapping Jx/Jy"
        if not runs:
            out.append({"x": float(xv), "y": None, "note": f"no pixels for {series} at t={xv:g}"})
            continue
        # Where dark & bright red separate around t=3..6, dark red (Jx) is slightly higher (smaller row)
        if series == "Jx_095":
            chosen = runs[0]
        elif series == "Jy_095":
            chosen = runs[-1]
        else:
            chosen = max(runs, key=lambda t: t[3])
        _, _, ctr, half = chosen
        slope = estimate_slope(used_mask, col, ctr, r_lo, r_hi)
        note = f"Fig3(c) {series}" + (f" ({note_extra})" if note_extra else "")
        out.append(make_reading(xv, pc, ctr, half, slope, xa, ya, extra_px=1.0, note=note))
    return out


# ---------------------------------------------------------------------------
# Figure 4 (e, f, g, h, i): Minimal oscillator trajectories and phase portraits
# ---------------------------------------------------------------------------

def read_fig4ef(image: np.ndarray, at: list[float], panel: str, color: str) -> list[dict]:
    arr = image.astype(np.float32)
    r, g, b = arr[:, :, 0], arr[:, :, 1], arr[:, :, 2]
    gray = arr.mean(axis=2)
    sat = arr.max(axis=2) - arr.min(axis=2)
    c_l, r_t, c_r, r_b = find_frame(image, thr=160.0, min_frac=0.45)
    ink = (gray < 170.0) & (sat < 35.0)

    bt = inward_ticks(ink, r_b[0], -1, c_l[1] + 5, c_r[0] - 5, is_x=True, min_len=6, max_len=40)
    lt = inward_ticks(ink, c_l[1], +1, r_t[1] + 5, r_b[0] - 5, is_x=False, min_len=6, max_len=40)
    if len(bt) < 2 or len(lt) < 2:
        raise RuntimeError(f"fig4({panel}) ticks failed: bt={bt}, lt={lt}")

    # Also detect the thin horizontal gray zero line across the plot interior!
    H_in = r_b[0] - r_t[1]
    zero_row_cands = []
    for row in range(r_t[1] + int(0.08 * H_in), r_b[0] - int(0.08 * H_in)):
        # Gray zero line has low saturation and gray < 210 across > 60% of box width
        frac = ((gray[row, c_l[1]+10:c_r[0]-10] < 225.0) & (sat[row, c_l[1]+10:c_r[0]-10] < 25.0)).mean()
        if frac > 0.55:
            zero_row_cands.append(row)
    z_gs = _group_indices(np.array(zero_row_cands, dtype=int), max_gap=3)

    if panel == "e":
        xa = LinearAxis([(bt[0][0], 0.0), (bt[-1][0], 60.0)])
        ymarks = [(lt[0][0], 0.1), (lt[-1][0], -0.4)]
        if z_gs:
            ymarks.append((z_gs[0][2], 0.0))
        ya = LinearAxis(ymarks)
    elif panel == "f":
        xa = LinearAxis([(bt[0][0], 0.0), (bt[-1][0], 140.0)])
        ymarks = [(lt[0][0], 0.2), (lt[-1][0], -0.2)]
        if z_gs:
            ymarks.append((z_gs[0][2], 0.0))
        ya = LinearAxis(ymarks)
    else:
        raise ValueError(panel)

    if color == "green":
        mask = (g - np.maximum(r, b) > 30.0) & (g > 80.0)
    else:
        mask = (b - np.maximum(r, g) > 30.0) & (b > 80.0)

    return read_branch_from_mask(
        mask, xa, ya, at, int(r_t[1] + 2), int(r_b[0] - 2),
        pick="thickest", search=4, extra_px=1.0, note_prefix=f"Fig4({panel}) {color}"
    )


def read_fig4ghi(image: np.ndarray, at: list[float], panel: str, branch: str) -> list[dict]:
    arr = image.astype(np.float32)
    r, g, b = arr[:, :, 0], arr[:, :, 1], arr[:, :, 2]
    gray = arr.mean(axis=2)
    sat = arr.max(axis=2) - arr.min(axis=2)
    c_l, r_t, c_r, r_b = find_frame(image, thr=170.0, min_frac=0.45)
    ink = (gray < 170.0) & (sat < 35.0)

    bt = inward_ticks(ink, r_b[0], -1, c_l[1] + 3, c_r[0] - 3, is_x=True, min_len=6, max_len=35)
    lt = inward_ticks(ink, c_l[1], +1, r_t[1] + 3, r_b[0] - 3, is_x=False, min_len=6, max_len=35)
    if len(bt) < 2:
        raise RuntimeError(f"fig4({panel}) bot ticks failed: {bt}")

    if panel == "g":
        xa = LinearAxis([(bt[0][0], -0.4), (bt[-1][0], 0.4)])
        if len(lt) < 2:
            raise RuntimeError(f"fig4g left ticks failed: {lt}")
        ya = LinearAxis([(lt[0][0], 0.4), (lt[-1][0], -0.4)])
    else:
        # Panels (h) and (i) share the exact [-0.4, 0.4] vertical scale with (g):
        # ticks are inset identically from top/bottom spines (0.74% of box height, matching (g)!)
        xa = LinearAxis([(bt[0][0], -0.3), (bt[-1][0], 0.3)])
        # Check if left ticks exist on shared spine; otherwise use exact proportional inset from spine
        if len(lt) >= 2 and (lt[-1][0] - lt[0][0]) > 0.9 * (r_b[2] - r_t[2]):
            ya = LinearAxis([(lt[0][0], 0.4), (lt[-1][0], -0.4)])
        else:
            # In panel (g), tick inset from spine centerline is 10.5 px on 1448.5 px = 0.00725 * H
            inset = 0.00725 * (r_b[2] - r_t[2])
            ya = LinearAxis([(r_t[2] + inset, 0.4), (r_b[2] - inset, -0.4)])

    # Pure bright red trajectory/limit cycle (excludes orange/green background streamlines!)
    pure_red = (r > 210.0) & (g < 50.0) & (b < 50.0)
    r_lo = int(r_t[1] + 8)
    r_hi = int(r_b[0] - 8)

    if panel == "g":
        # Lower-left red trajectory from (0,0) to stable node (-0.236, -0.311)
        r_mid = int(round(ya.to_pixel(0.02)))
        return read_branch_from_mask(
            pure_red, xa, ya, at, r_mid, r_hi, pick="thickest", search=4,
            extra_px=1.0, note_prefix="Fig4(g) red trajectory"
        )
    else:
        # Outer limit cycle: top = highest red crossing, bottom = lowest red crossing!
        pick = "top" if branch == "upper" else "bottom"
        return read_branch_from_mask(
            pure_red, xa, ya, at, r_lo, r_hi, pick=pick, search=3,
            min_run_px=3, extra_px=1.5, note_prefix=f"Fig4({panel}) limit cycle {branch}"
        )


# ---------------------------------------------------------------------------
# Figure 5 (a, b): Exceptional-point red boundaries in (chi, alpha) plane
# ---------------------------------------------------------------------------

def read_fig5(image: np.ndarray, at: list[float], panel: str) -> list[dict]:
    arr = image.astype(np.float32)
    r, g, b = arr[:, :, 0], arr[:, :, 1], arr[:, :, 2]
    gray = arr.mean(axis=2)
    sat = arr.max(axis=2) - arr.min(axis=2)
    c_l, r_t, c_r, r_b = find_frame(image, thr=160.0, min_frac=0.45)
    ink = (gray < 160.0) & (sat < 35.0)

    bt = inward_ticks(ink, r_b[0], -1, c_l[1] + 15, c_r[0] - 15, is_x=True, min_len=18, max_len=50)
    lt = inward_ticks(ink, c_l[1], +1, r_t[1] + 15, r_b[0] - 15, is_x=False, min_len=18, max_len=50)
    if len(bt) < 2 or len(lt) < 2:
        raise RuntimeError(f"fig5({panel}) ticks failed: bt={bt}, lt={lt}")

    xa = LinearAxis([(bt[0][0], -0.25), (bt[-1][0], 0.25)])
    ymax = 0.33 if panel == "a" else 0.25
    ya = LinearAxis([(lt[0][0], ymax), (lt[-1][0], -ymax)])

    red = (r > 200.0) & (g < 75.0) & (b < 75.0)
    r_lo = int(r_t[1] + 10)
    r_mid = int(round(ya.to_pixel(-0.005)))
    return read_branch_from_mask(
        red, xa, ya, at, r_lo, r_mid, pick="top", search=3,
        extra_px=1.5, note_prefix=f"Fig5({panel}) upper EP boundary"
    )


# ---------------------------------------------------------------------------
# Figure 6 (b, c, d): Hopf circle C_r, hyperbola C_i, blue arm, and inner boundaries
# ---------------------------------------------------------------------------

def read_fig6b(image: np.ndarray, at: list[float], branch: str) -> list[dict]:
    arr = image.astype(np.float32)
    r, g, b = arr[:, :, 0], arr[:, :, 1], arr[:, :, 2]
    gray = arr.mean(axis=2)
    sat = arr.max(axis=2) - arr.min(axis=2)
    c_l, r_t, c_r, r_b = find_frame(image, thr=160.0, min_frac=0.45)
    ink = (gray < 160.0) & (sat < 35.0)

    bt = inward_ticks(ink, r_b[0], -1, c_l[1] + 15, c_r[0] - 15, is_x=True, min_len=15, max_len=45)
    lt = inward_ticks(ink, c_l[1], +1, r_t[1] + 15, r_b[0] - 15, is_x=False, min_len=15, max_len=45)
    if len(bt) != 3 or len(lt) != 3:
        raise RuntimeError(f"fig6b expected 3 x and 3 y ticks, got bt={bt}, lt={lt}")
    xa = LinearAxis([(bt[0][0], 0.0), (bt[1][0], 0.5), (bt[2][0], 1.0)])
    ya = LinearAxis([(lt[0][0], 1.0), (lt[1][0], 0.5), (lt[2][0], 0.0)])

    red = (r > 155.0) & (r - g > 45.0) & (r - b > 45.0)
    r_mid = int(round(ya.to_pixel(0.30)))
    if branch == "upper":
        r_lo = int(round(ya.to_pixel(0.48)))
        r_hi = r_mid
    else:
        r_lo = r_mid
        r_hi = int(round(ya.to_pixel(0.12)))
    pick = "top" if branch == "upper" else "bottom"
    return read_branch_from_mask(
        red, xa, ya, at, r_lo, r_hi, pick=pick, search=3,
        extra_px=1.0, note_prefix=f"Fig6(b) red C_r circle {branch}"
    )


def _calibrate_fig6c(image: np.ndarray):
    arr = image.astype(np.float32)
    gray = arr.mean(axis=2)
    c_l, r_t, c_r, r_b = find_frame(image, thr=140.0, min_frac=0.45)
    # Use dark luminance < 120 (without saturation cut) so the 0.5 tick next to the blue arm is detected cleanly
    dark = gray < 120.0
    tt = inward_ticks(dark, r_t[1], +1, c_l[1] + 20, c_r[0] - 20, is_x=True, min_len=28, max_len=65)
    bt = inward_ticks(dark, r_b[0], -1, c_l[1] + 20, c_r[0] - 20, is_x=True, min_len=28, max_len=65)
    if len(tt) != 3 or len(bt) != 3:
        raise RuntimeError(f"fig6c expected 3 x ticks on top/bot, got tt={tt}, bt={bt}")
    x3 = [0.5 * (t[0] + b[0]) for t, b in zip(tt, bt)]
    xa = LinearAxis([(x3[0], 0.0), (x3[1], 0.5), (x3[2], 1.0)])
    # Right spine ticks: 7 major ticks at [1.0, 0.8, 0.6, 0.4, 0.2, 0.0, -0.2] (~292.4 px apart)
    rt = inward_ticks(dark, c_r[0], -1, r_t[1] + 15, r_b[0] - 15, is_x=False, min_len=35, max_len=65)
    step_est = (r_b[2] - r_t[2]) / 7.06
    ys = [p for p, _ in rt]
    y7 = [ys[0]]
    for k in range(1, 7):
        target = y7[0] + k * step_est
        cand = min(ys, key=lambda p: abs(p - target))
        if abs(cand - target) < 25.0:
            y7.append(cand)
    if len(y7) != 7:
        raise RuntimeError(f"fig6c expected 7 y ticks, got {y7} from {rt}")
    ya = LinearAxis(list(zip(y7, [1.0, 0.8, 0.6, 0.4, 0.2, 0.0, -0.2])))
    return xa, ya, (c_l, r_t, c_r, r_b)


def read_fig6c(image: np.ndarray, at: list[float], branch: str) -> list[dict]:
    xa, ya, (c_l, r_t, c_r, r_b) = _calibrate_fig6c(image)
    arr = image.astype(np.float32)
    r, g, b = arr[:, :, 0], arr[:, :, 1], arr[:, :, 2]
    r_lo = int(r_t[1] + 6)
    r_hi = int(r_b[0] - 6)
    r_mid = int(round(ya.to_pixel(0.30)))

    if branch == "ci_top":
        orange = (r > 180.0) & (r - b > 80.0) & (g < 180.0)
        return read_branch_from_mask(
            orange, xa, ya, at, r_lo, r_mid, pick="top", search=3,
            extra_px=1.0, note_prefix="Fig6(c) orange C_i top"
        )
    elif branch == "ci_bot":
        orange = (r > 180.0) & (r - b > 80.0) & (g < 180.0)
        return read_branch_from_mask(
            orange, xa, ya, at, r_mid, r_hi, pick="bottom", search=3,
            extra_px=1.0, note_prefix="Fig6(c) orange C_i bottom"
        )
    elif branch == "ci_right_upper":
        orange = (r > 180.0) & (r - b > 80.0) & (g < 180.0)
        return read_branch_from_mask(
            orange, xa, ya, at, r_lo, r_mid, pick="bottom", search=3,
            extra_px=1.0, note_prefix="Fig6(c) orange C_i right-upper"
        )
    elif branch == "blue_arm_right_upper":
        # Right blue spinodal arm: inner solid blue boundary in phi2 in [0.30, 0.425]
        # (strictly below the outer binodal at phi2 ~0.44-0.49 for phi1 >= 0.65).
        blue = (b - r > 30.0) & (b > 110.0) & (r < 130.0) & (g < 190.0)
        r_arm_top = int(round(ya.to_pixel(0.425)))
        r_arm_mid = int(round(ya.to_pixel(0.30)))
        return read_branch_from_mask(
            blue, xa, ya, at, r_arm_top, r_arm_mid, pick="top", search=4,
            min_run_px=1, extra_px=1.0, note_prefix="Fig6(c) right blue spinodal arm upper"
        )
    raise ValueError(f"Unknown fig6c branch {branch}")


def read_fig6d(image: np.ndarray, at: list[float], curve: str, branch: str) -> list[dict]:
    arr = image.astype(np.float32)
    r, g, b = arr[:, :, 0], arr[:, :, 1], arr[:, :, 2]
    gray = arr.mean(axis=2)
    sat = arr.max(axis=2) - arr.min(axis=2)
    c_l, r_t, c_r, r_b = find_frame(image, thr=140.0, min_frac=0.45)
    ink = (gray < 140.0) & (sat < 35.0)

    tt = inward_ticks(ink, r_t[1], +1, c_l[1] + 20, c_r[0] - 20, is_x=True, min_len=28, max_len=65)
    lt = inward_ticks(ink, c_l[1], +1, r_t[1] + 20, r_b[0] - 20, is_x=False, min_len=28, max_len=65)
    if len(tt) != 4 or len(lt) != 3:
        raise RuntimeError(f"fig6d ticks failed: tt={tt}, lt={lt}")
    xa = LinearAxis([(tt[0][0], 0.2), (tt[1][0], 0.3), (tt[2][0], 0.4), (tt[3][0], 0.5)])
    ya = LinearAxis([(r_t[2], 0.5), (lt[0][0], 0.4), (lt[1][0], 0.3), (lt[2][0], 0.2), (r_b[2], 0.1)])

    r_lo = int(r_t[1] + 15)
    r_hi = int(r_b[0] - 15)
    r_mid = int(round(ya.to_pixel(0.30)))
    pick = "top" if branch == "upper" else "bottom"

    if curve == "black_outer":
        mask = (r < 60.0) & (g < 60.0) & (b < 60.0)
        # Mask out the colorbar on the far right (phi1 > 0.505)
        c_bar = int(round(xa.to_pixel(0.508)))
        mask[:, c_bar:] = False
        r_a, r_b_lim = (r_lo, r_mid) if branch == "upper" else (r_mid, r_hi)
        return read_branch_from_mask(
            mask, xa, ya, at, r_a, r_b_lim, pick=pick, search=3,
            min_run_px=3, extra_px=1.0, note_prefix=f"Fig6(d) black outer {branch}"
        )
    elif curve == "gray_inner":
        # Medium gray curve (R, G, B in 105..195, very low saturation < 18)
        mask = (gray > 105.0) & (gray < 195.0) & (sat < 18.0)
        # Restrict strictly to inner gray loop region (phi1 in [0.25, 0.45], phi2 in [0.18, 0.42])
        r_top_g = int(round(ya.to_pixel(0.42)))
        r_bot_g = int(round(ya.to_pixel(0.18)))
        r_a, r_b_lim = (r_top_g, r_mid) if branch == "upper" else (r_mid, r_bot_g)
        return read_branch_from_mask(
            mask, xa, ya, at, r_a, r_b_lim, pick=pick, search=5,
            min_run_px=2, extra_px=1.5, note_prefix=f"Fig6(d) gray inner {branch}"
        )
    raise ValueError(f"Unknown fig6d curve={curve}")
