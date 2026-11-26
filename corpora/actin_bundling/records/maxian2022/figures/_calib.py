"""Shared image-based axis calibration and figure reading helpers for maxian2022.

Every axis calibration and data reading is derived strictly from the passed RGB
`image` array (no PDF re-rendering inside `read(image, at)` and no hardcoded
output values). Pure vectorized NumPy operations keep every reading < 0.1 s.
"""

from __future__ import annotations

import math
import numpy as np

from boxgraph import digitize as D


# ---------------------------------------------------------------------------
# Standard RGB colors from Maxian et al. (2022) MATLAB vector plots
# ---------------------------------------------------------------------------
MATLAB_BLUE = (0, 114, 189)
MATLAB_ORANGE = (217, 83, 25)
MATLAB_YELLOW = (237, 177, 32)
MATLAB_PURPLE = (126, 47, 142)

FIG7_MAROON = (162, 20, 47)      # B, tau_f = 2 tau_c
FIG7_PINK = (255, 153, 204)      # NB, tau_f = 2 tau_c
FIG7_DGREEN = (0, 128, 0)        # B, tau_f = tau_c
FIG7_BGREEN = (20, 255, 0)       # NB, tau_f = tau_c
FIG7_BLACK = (0, 0, 0)           # B, tau_f = tau_c / 2
FIG7_GRAY = (148, 148, 148)      # NB, tau_f = tau_c / 2

FIGS1_BAR_BLUE = (100, 172, 220)

FIGS7_DBLUE = (0, 114, 189)      # lm=0.2, L=1
FIGS7_LBLUE = (77, 190, 238)     # lm=0.2, L=0.5
FIGS7_RED = (255, 0, 0)          # lm=0.4, L=1
FIGS7_MAGENTA = (255, 0, 255)    # lm=0.4, L=0.5
FIGS7_BLACK = (0, 0, 0)          # lm=0.8, L=1
FIGS7_LGRAY = (204, 204, 204)    # lm=0.8, L=0.5


# ---------------------------------------------------------------------------
# Frame & tick detection (fast vectorized uint8 / strip-sliced operations)
# ---------------------------------------------------------------------------

def _longest_run_1d(mask_1d: np.ndarray, bridge: int = 2) -> tuple[int, int, int]:
    """Fast vectorized longest True run in a 1D boolean array."""
    if not mask_1d.any():
        return 0, 0, 0
    m = mask_1d.copy()
    if bridge >= 1 and len(m) > 2:
        m[1:-1] |= (m[:-2] & m[2:])
    if bridge >= 2 and len(m) > 3:
        m[1:-2] |= (m[:-3] & m[3:])
        m[2:-1] |= (m[:-3] & m[3:])
    d = np.diff(np.concatenate(([0], m.view(np.uint8), [0])))
    s = np.nonzero(d == 1)[0]
    e = np.nonzero(d == -1)[0]
    if len(s) == 0:
        return 0, 0, 0
    lens = e - s
    idx = int(np.argmax(lens))
    return int(lens[idx]), int(s[idx]), int(e[idx])


def find_outer_frame(
    img: np.ndarray,
    *,
    min_w_frac: float = 0.55,
    min_h_frac: float = 0.55,
    dark_thr: int = 215,
) -> tuple[float, float, float, float]:
    """Locate the outermost rectangular MATLAB plot box (cL, rT, cR, rB)."""
    # Fast uint8 non-white test (no full-image float32 allocation!)
    nw = img.min(axis=2) < dark_thr
    H, W = nw.shape

    row_sums = nw.sum(axis=1)
    h_cands: list[tuple[int, int, int]] = []
    for r in np.nonzero(row_sums >= int(min_w_frac * W))[0]:
        ln, s, e = _longest_run_1d(nw[r], bridge=2)
        if ln >= min_w_frac * W:
            h_cands.append((int(r), s, e))
    if not h_cands:
        raise RuntimeError("No horizontal frame spines found")

    max_w = max(e - s for _, s, e in h_cands)
    wide_h = [x for x in h_cands if (x[2] - x[1]) >= 0.90 * max_w]
    r_min = wide_h[0][0]
    r_max = wide_h[-1][0]
    top_rows = [r for r, _, _ in wide_h if r <= r_min + 6]
    bot_rows = [r for r, _, _ in wide_h if r >= r_max - 6]
    rT = float(np.mean(top_rows)) + 0.5
    rB = float(np.mean(bot_rows)) + 0.5

    r0i = max(0, int(rT) - 2)
    r1i = min(H, int(rB) + 3)
    sub_h = r1i - r0i
    col_sums = nw[r0i:r1i, :].sum(axis=0)
    v_cands: list[int] = []
    for c in np.nonzero(col_sums >= int(0.88 * sub_h))[0]:
        ln, _, _ = _longest_run_1d(nw[r0i:r1i, c], bridge=2)
        if ln >= 0.88 * sub_h:
            v_cands.append(int(c))
    if not v_cands:
        raise RuntimeError("No vertical frame spines found")
    c_min = v_cands[0]
    c_max = v_cands[-1]
    left_cols = [c for c in v_cands if c <= c_min + 6]
    right_cols = [c for c in v_cands if c >= c_max - 6]
    cL = float(np.mean(left_cols)) + 0.5
    cR = float(np.mean(right_cols)) + 0.5
    return cL, rT, cR, rB


def find_inward_ticks(
    img: np.ndarray,
    frame: tuple[float, float, float, float],
    side: str,
    *,
    min_len: int = 8,
    max_len: int = 45,
    start_off: int = 4,
    gray_thr: int = 145,
    sat_thr: int = 32,
) -> list[tuple[float, int]]:
    """Find inward neutral-dark tick marks along one spine ('bottom','top','left','right').

    Slices only the thin 45-px border strip before thresholding (~1 ms total).
    Returns sorted list of (center_px, tick_length_px).
    """
    cL, rT, cR, rB = frame
    H, W = img.shape[:2]

    if side in ("bottom", "top"):
        c0 = max(0, int(round(cL)) + 5)
        c1 = min(W, int(round(cR)) - 5)
        if side == "bottom":
            r_start = int(round(rB)) - start_off
            r_end = max(0, r_start - max_len)
            sub = img[r_end:r_start + 1, c0:c1][::-1, :, :]  # axis 0 goes inward
        else:
            r_start = int(round(rT)) + start_off
            r_end = min(H, r_start + max_len)
            sub = img[r_start:r_end, c0:c1, :]
        offset = c0
    else:
        r0 = max(0, int(round(rT)) + 5)
        r1 = min(H, int(round(rB)) - 5)
        if side == "left":
            c_start = int(round(cL)) + start_off
            c_end = min(W, c_start + max_len)
            sub = img[r0:r1, c_start:c_end, :].transpose(1, 0, 2)  # axis 0 goes inward
        else:
            c_start = int(round(cR)) - start_off
            c_end = max(0, c_start - max_len)
            sub = img[r0:r1, c_end:c_start + 1, :][:, ::-1, :].transpose(1, 0, 2)
        offset = r0

    if sub.shape[0] < min_len or sub.shape[1] == 0:
        return []

    # Neutral dark test ONLY on the tiny sliced strip
    mx = sub.max(axis=2).astype(np.int16)
    mn = sub.min(axis=2).astype(np.int16)
    nd = ((mx + mn) < (2 * gray_thr)) & ((mx - mn) < sat_thr)

    # Vectorized run length starting within first 2 inward pixels
    # Pad one False row at the end so argmax(~nd) always finds the end of the run
    L, N = nd.shape
    padded = np.vstack([nd, np.zeros((1, N), dtype=bool)])
    # Start index k in {0, 1, 2} where True first appears in first 2 rows
    s0 = padded[0]
    s1 = (~s0) & padded[1]
    active = s0 | s1
    if not active.any():
        return []
    hits: list[tuple[int, int]] = []
    for j in np.nonzero(active)[0]:
        k = 0 if s0[j] else 1
        col = padded[k:, j]
        ln = int(np.argmax(~col))
        if min_len <= ln <= max_len - 2:
            hits.append((offset + int(j), ln))

    grps: list[list[tuple[int, int]]] = []
    for p, ln in hits:
        if not grps or p > grps[-1][-1][0] + 4:
            grps.append([(p, ln)])
        else:
            grps[-1].append((p, ln))
    out: list[tuple[float, int]] = []
    for g in grps:
        if len(g) <= 9:
            center = 0.5 * (g[0][0] + g[-1][0]) + 0.5
            out.append((center, max(x[1] for x in g)))
    return out


def calibrate_matlab_axes(
    img: np.ndarray,
    *,
    x_bounds: tuple[float, float] | None = None,
    y_bounds: tuple[float, float] | None = None,
    x_scale: str = "linear",
    y_scale: str = "linear",
    x_interior_vals: list[float] | None = None,
    y_interior_vals: list[float] | None = None,
    x_tick_side: str = "top",
    y_tick_side: str = "left",
    x_all_ticks: list[float] | None = None,
    y_all_ticks: list[float] | None = None,
    major_only_x: bool = True,
    major_only_y: bool = True,
    max_res_px: float = 2.0,
) -> tuple[D.Axis, D.Axis, tuple[float, float, float, float]]:
    """Calibrate x and y axes from the image frame and inward tick marks.

    When multiple interior ticks are present, fits directly from detected ticks
    and only adds frame borders when they lie on the same scale (within 2.5 px),
    automatically handling MATLAB plots whose y_min/y_max has small auto-padding
    (such as Fig. S5 y_min ~ -0.015 or Fig. 3 alignment y_min ~ 0.698).
    """
    frame = find_outer_frame(img)
    cL, rT, cR, rB = frame

    def _fit_one(
        scale: str,
        bounds: tuple[float, float] | None,
        interior_vals: list[float] | None,
        all_ticks: list[float] | None,
        side: str,
        p_lo: float,
        p_hi: float,
        major_only: bool,
        invert: bool,
    ) -> D.Axis:
        opp = {"top": "bottom", "bottom": "top", "left": "right", "right": "left"}[side]
        raw1 = find_inward_ticks(img, frame, side)
        raw2 = find_inward_ticks(img, frame, opp)

        def _filter_major(raw):
            if not raw:
                return []
            if not major_only:
                return [p for p, _ in raw]
            mx = max(ln for _, ln in raw)
            return [p for p, ln in raw if ln >= 0.72 * mx]

        cand1 = _filter_major(raw1)
        cand2 = _filter_major(raw2)

        if all_ticks is not None:
            cand = cand1 if len(cand1) == len(all_ticks) else (cand2 if len(cand2) == len(all_ticks) else cand1)
            if len(cand) != len(all_ticks):
                raise RuntimeError(
                    f"Expected {len(all_ticks)} ticks on {side}/{opp}, found {len(cand1)}/{len(cand2)}: {cand1}"
                )
            vals = list(reversed(all_ticks)) if invert else list(all_ticks)
            ax = D.Axis(scale, list(zip(cand, vals)))
            if ax.residual_px > max_res_px:
                raise RuntimeError(f"Axis residual {ax.residual_px:.2f} px > {max_res_px}")
            return ax

        assert bounds is not None
        v0, v1 = bounds
        v_at_lo = v1 if invert else v0
        v_at_hi = v0 if invert else v1
        border_marks = [(p_lo, v_at_lo), (p_hi, v_at_hi)]

        if interior_vals:
            # Preliminary linear/log map from borders (good within ~20 px even with small MATLAB padding)
            prelim = D.Axis(scale, border_marks)
            pool1 = [p for p, _ in raw1]
            pool2 = [p for p, _ in raw2]

            # Also check if v0 or v1 themselves appear as inward ticks slightly inside the border!
            all_expected = list(interior_vals)
            for vb in (v0, v1):
                if vb not in all_expected:
                    all_expected.append(vb)

            matched_interior: list[tuple[float, float]] = []
            matched_boundary_tick: dict[float, float] = {}
            for v in all_expected:
                p_exp = prelim.to_pixel(v)
                tol_px = 20.0 if v in (v0, v1) else 15.0
                best = None
                for pool in (pool1, pool2):
                    if pool:
                        near = min(pool, key=lambda p: abs(p - p_exp))
                        if abs(near - p_exp) <= tol_px:
                            best = near
                            break
                if best is not None:
                    if v in (v0, v1):
                        matched_boundary_tick[v] = best
                    else:
                        matched_interior.append((best, v))

            # If we matched >= 2 interior ticks cleanly, trust the ticks first!
            if len(matched_interior) >= 2:
                tick_marks = list(matched_interior)
                # Check if boundary ticks exist and are consistent with interior ticks
                ax_int = D.Axis(scale, tick_marks)
                if ax_int.residual_px <= max_res_px:
                    for vb, p_border in ((v_at_lo, p_lo), (v_at_hi, p_hi)):
                        p_pred = ax_int.to_pixel(vb)
                        if vb in matched_boundary_tick and abs(matched_boundary_tick[vb] - p_pred) <= 2.5:
                            tick_marks.append((matched_boundary_tick[vb], vb))
                        elif abs(p_border - p_pred) <= 2.5:
                            tick_marks.append((p_border, vb))
                    ax = D.Axis(scale, tick_marks)
                    if ax.residual_px <= max_res_px:
                        return ax

            # Fallback: combine borders with closely matched interior ticks
            marks = list(border_marks)
            for p, v in matched_interior:
                if abs(p - prelim.to_pixel(v)) <= 5.0:
                    marks.append((p, v))
            ax = D.Axis(scale, marks)
            if ax.residual_px > max_res_px:
                raise RuntimeError(f"Axis residual {ax.residual_px:.2f} px > {max_res_px}")
            return ax

        ax = D.Axis(scale, border_marks)
        return ax

    xa = _fit_one(x_scale, x_bounds, x_interior_vals, x_all_ticks, x_tick_side, cL, cR, major_only_x, invert=False)
    ya = _fit_one(y_scale, y_bounds, y_interior_vals, y_all_ticks, y_tick_side, rT, rB, major_only_y, invert=True)
    return xa, ya, frame


def calibrate_log_y_from_decades(
    img: np.ndarray,
    frame: tuple[float, float, float, float],
    major_decade_vals: list[float],
    *,
    side: str = "left",
    max_res_px: float = 1.5,
) -> D.Axis:
    """Calibrate a log10 y-axis whose frame borders are NOT powers of 10 (e.g., Fig. 2 bottom-right)."""
    raw = find_inward_ticks(img, frame, side, min_len=6)
    if not raw:
        raw = find_inward_ticks(img, frame, "right", min_len=6)
    mx = max(ln for _, ln in raw)
    majors = [p for p, ln in raw if ln >= 0.72 * mx]
    if len(majors) != len(major_decade_vals):
        raise RuntimeError(f"Expected {len(major_decade_vals)} major log y-ticks, found {len(majors)}: {majors}")
    prelim = D.Axis("log10", list(zip(majors, major_decade_vals)))
    marks = list(zip(majors, major_decade_vals))
    min_d = int(math.floor(math.log10(min(major_decade_vals)))) - 1
    max_d = int(math.ceil(math.log10(max(major_decade_vals)))) + 1
    for p, _ in raw:
        for k in range(min_d, max_d + 1):
            for m in range(2, 10):
                vm = m * (10.0 ** k)
                if abs(prelim.to_pixel(vm) - p) <= 2.5:
                    marks.append((p, vm))
                    break
    ax = D.Axis("log10", marks)
    if ax.residual_px > max_res_px:
        raise RuntimeError(f"Log y-axis residual {ax.residual_px:.2f} px > {max_res_px}")
    return ax


# ---------------------------------------------------------------------------
# Color masking & exclusion regions (fast uint8/int16 inside frame only)
# ---------------------------------------------------------------------------

def rgb_mask(
    img: np.ndarray,
    rgb: tuple[int, int, int],
    *,
    tol: float = 48.0,
    frame: tuple[float, float, float, float] | None = None,
    inset: int = 4,
    exclude_fracs: list[tuple[float, float, float, float]] = (),
    black_mode: bool = False,
) -> np.ndarray:
    """Boolean mask for pixels matching `rgb` inside `frame` minus `exclude_fracs`."""
    H, W = img.shape[:2]
    out = np.zeros((H, W), dtype=bool)
    if frame is not None:
        cL, rT, cR, rB = frame
        r0 = max(0, int(math.ceil(rT)) + inset)
        r1 = min(H, int(math.floor(rB)) - inset)
        c0 = max(0, int(math.ceil(cL)) + inset)
        c1 = min(W, int(math.floor(cR)) - inset)
    else:
        r0, r1, c0, c1 = 0, H, 0, W

    sub = img[r0:r1, c0:c1]
    if black_mode:
        m = sub.max(axis=2) < 65
    else:
        dr = sub[:, :, 0].astype(np.int32) - int(rgb[0])
        dg = sub[:, :, 1].astype(np.int32) - int(rgb[1])
        db = sub[:, :, 2].astype(np.int32) - int(rgb[2])
        m = (dr * dr + dg * dg + db * db) <= int(tol * tol)
    out[r0:r1, c0:c1] = m

    if frame is not None:
        cL, rT, cR, rB = frame
        for fx0, fy0, fx1, fy1 in exclude_fracs:
            ec0 = int(round(cL + fx0 * (cR - cL)))
            ec1 = int(round(cL + fx1 * (cR - cL)))
            er0 = int(round(rT + fy0 * (rB - rT)))
            er1 = int(round(rT + fy1 * (rB - rT)))
            out[max(0, er0):min(H, er1), max(0, ec0):min(W, ec1)] = False
    return out


# ---------------------------------------------------------------------------
# Reading Primitive 1: MATLAB errorbar time series (mean + SEM bars)
# ---------------------------------------------------------------------------

def _vertical_runs_bridged(
    col_mask: np.ndarray,
    col_other_ink: np.ndarray,
    max_ink_gap: int = 35,
    max_white_gap: int = 2,
) -> list[tuple[int, int, int]]:
    """Return vertical runs (r_start, r_end, true_pixel_count) in one column.

    Bridges tiny gaps (<= max_white_gap) AND bridges gaps up to `max_ink_gap`
    ONLY when another curve's ink occupies the gap.
    """
    ys = np.nonzero(col_mask)[0]
    if len(ys) == 0:
        return []
    runs: list[tuple[int, int, int]] = []
    rs = int(ys[0])
    re = rs
    cnt = 1
    for y in ys[1:]:
        yi = int(y)
        gap = yi - re - 1
        bridge_ok = False
        if gap <= max_white_gap:
            bridge_ok = True
        elif gap <= max_ink_gap:
            if col_other_ink[re + 1:yi].mean() >= 0.70:
                bridge_ok = True
        if bridge_ok:
            re = yi
            cnt += 1
        else:
            runs.append((rs, re, cnt))
            rs = re = yi
            cnt = 1
    runs.append((rs, re, cnt))
    return runs


def read_errorbar_series(
    img: np.ndarray,
    xa: D.Axis,
    ya: D.Axis,
    frame: tuple[float, float, float, float],
    rgb: tuple[int, int, int],
    at: list[float],
    *,
    tol: float = 48.0,
    exclude_fracs: list[tuple[float, float, float, float]] = (),
    black_mode: bool = False,
    stem_search_px: int = 5,
    cap_half_px: float = 3.5,
    expected_y_band: tuple[float, float] | None = None,
    y_guide_fn=None,
    guide_half_px: float | None = None,
    smooth_fallback: bool = True,
    max_ink_gap: int = 35,
    min_run_px: int = 3,
    fallback_search_px: int = 12,
    use_center_when_short: bool = True,
    blank_ticks: bool = False,
) -> list[dict]:
    """Read a MATLAB `errorbar(t, y, sem)` series at the x positions in `at`."""
    cL, rT, cR, rB = frame
    mask = rgb_mask(img, rgb, tol=tol, frame=frame, inset=3, exclude_fracs=exclude_fracs, black_mode=black_mode)
    if blank_ticks or black_mode:
        rTi, rBi, cLi, cRi = int(round(rT)), int(round(rB)), int(round(cL)), int(round(cR))
        mask[rTi:rTi + 28, :] = False
        mask[max(0, rBi - 28):rBi + 1, :] = False
        mask[:, cLi:cLi + 28] = False
        mask[:, max(0, cRi - 28):cRi + 1] = False

    other_ink = (img.min(axis=2) < 225) & (~mask)

    r_lo_lim = int(math.ceil(rT)) + 3
    r_hi_lim = int(math.floor(rB)) - 3
    if expected_y_band is not None:
        p1 = ya.to_pixel(expected_y_band[0])
        p2 = ya.to_pixel(expected_y_band[1])
        r_lo_lim = max(r_lo_lim, int(min(p1, p2)))
        r_hi_lim = min(r_hi_lim, int(max(p1, p2)))

    cal_unc_px = max(ya.residual_px, 0.5) + 2.0
    out: list[dict] = []

    for x in at:
        c_nom = xa.to_pixel(x)
        c_int = int(round(c_nom))
        r_min_use, r_max_use = r_lo_lim, r_hi_lim
        if y_guide_fn is not None and guide_half_px is not None:
            r_g = ya.to_pixel(y_guide_fn(x))
            r_min_use = max(r_lo_lim, int(round(r_g - guide_half_px)))
            r_max_use = min(r_hi_lim, int(round(r_g + guide_half_px)))

        def _best_in_window(half_w: int) -> tuple[int, int, int, int] | None:
            c0 = max(int(math.ceil(cL)) + 3, c_int - half_w)
            c1 = min(int(math.floor(cR)) - 3, c_int + half_w)
            best = None
            for c in range(c0, c1 + 1):
                runs = _vertical_runs_bridged(
                    mask[r_min_use:r_max_use, c],
                    other_ink[r_min_use:r_max_use, c],
                    max_ink_gap=max_ink_gap,
                )
                for rs, re, cnt in runs:
                    if cnt < min_run_px:
                        continue
                    rs_a = rs + r_min_use
                    re_a = re + r_min_use
                    span = re_a - rs_a + 1
                    dist_pen = abs(c - c_nom) * 1.5
                    score = cnt + 0.6 * span - dist_pen
                    if best is None or score > best[0]:
                        best = (score, c, rs_a, re_a, cnt)
            if best is None:
                return None
            return best[1], best[2], best[3], best[4]

        hit = _best_in_window(stem_search_px)
        note_extra = ""
        if hit is None and smooth_fallback:
            hit = _best_in_window(fallback_search_px)
            if hit is not None:
                note_extra = " (nearby col)"

        if hit is None:
            out.append({"x": float(x), "y": None, "found": False, "note": "no matching color pixels near x"})
            continue

        c_best, rs_a, re_a, cnt = hit
        span = re_a - rs_a + 1
        clipped_top = rs_a <= int(round(rT)) + 6
        clipped_bot = re_a >= int(round(rB)) - 6

        if span <= 16 and use_center_when_short:
            r_mid = 0.5 * (rs_a + re_a) + 0.5
            r_hi_px = float(rs_a) + 0.5
            r_lo_px = float(re_a) + 0.5
            y_val = ya.to_value(r_mid)
            bar_hi = ya.to_value(r_hi_px)
            bar_lo = ya.to_value(r_lo_px)
        else:
            r_hi_px = float(rs_a) + (0.5 if clipped_top else cap_half_px)
            r_lo_px = float(re_a) + (0.5 if clipped_bot else (1.0 - cap_half_px))
            if r_lo_px < r_hi_px:
                r_hi_px = r_lo_px = 0.5 * (rs_a + re_a) + 0.5
            bar_hi = ya.to_value(r_hi_px)
            bar_lo = ya.to_value(r_lo_px)
            if clipped_bot and not clipped_top:
                r_mid_cands = []
                for dc in (-9, -8, -7, 7, 8, 9):
                    cc = c_best + dc
                    if int(cL) + 4 <= cc <= int(cR) - 4:
                        ys_c = np.nonzero(mask[rs_a:re_a + 1, cc])[0]
                        if len(ys_c) >= 3:
                            r_mid_cands.append(float(np.median(ys_c)) + rs_a + 0.5)
                if r_mid_cands:
                    r_mid = float(np.median(r_mid_cands))
                    y_val = ya.to_value(r_mid)
                    note_extra += " (bottom bar clipped)"
                else:
                    y_val = 0.5 * (bar_lo + bar_hi)
            elif clipped_top and not clipped_bot:
                r_mid_cands = []
                for dc in (-9, -8, -7, 7, 8, 9):
                    cc = c_best + dc
                    if int(cL) + 4 <= cc <= int(cR) - 4:
                        ys_c = np.nonzero(mask[rs_a:re_a + 1, cc])[0]
                        if len(ys_c) >= 3:
                            r_mid_cands.append(float(np.median(ys_c)) + rs_a + 0.5)
                if r_mid_cands:
                    r_mid = float(np.median(r_mid_cands))
                    y_val = ya.to_value(r_mid)
                    note_extra += " (top bar clipped)"
                else:
                    y_val = 0.5 * (bar_lo + bar_hi)
            else:
                y_val = 0.5 * (bar_lo + bar_hi)

        b_lo, b_hi = min(bar_lo, bar_hi), max(bar_lo, bar_hi)
        y_val = max(b_lo, min(b_hi, y_val))
        row_y = ya.to_pixel(y_val)

        y_u1 = ya.to_value(row_y + cal_unc_px)
        y_u2 = ya.to_value(row_y - cal_unc_px)
        lo_val = min(y_u1, y_u2, y_val)
        hi_val = max(y_u1, y_u2, y_val)

        bh_out = max(b_hi, ya.to_value(r_hi_px - cal_unc_px))
        bl_out = min(b_lo, ya.to_value(r_lo_px + cal_unc_px))
        if ya.scale == "linear" and b_lo >= 0.0 and bl_out < 0.0 and y_val >= 0.0:
            if abs(ya.to_value(rB)) < 0.05 * abs(ya.to_value(rT) - ya.to_value(rB)):
                bl_out = max(0.0, bl_out)

        out.append({
            "x": float(x),
            "y": float(y_val),
            "lo": float(lo_val),
            "hi": float(hi_val),
            "bar_lo": float(b_lo),
            "bar_hi": float(b_hi),
            "bar_lo_outer": float(bl_out),
            "bar_hi_outer": float(bh_out),
            "px": [float(c_best) + 0.5, float(row_y)],
            "px_lo": float(r_hi_px),
            "px_hi": float(r_lo_px),
            "note": f"stem {span}px (px_cnt={cnt}){note_extra}",
        })
    return out


# ---------------------------------------------------------------------------
# Reading Primitive 2: Dashed/dotted black fit curves
# ---------------------------------------------------------------------------

def read_black_fit_curve(
    img: np.ndarray,
    xa: D.Axis,
    ya: D.Axis,
    frame: tuple[float, float, float, float],
    at: list[float],
    *,
    branch: str = "only",
    branch_split_x: float | None = None,
    exclude_fracs: list[tuple[float, float, float, float]] = (),
    half_win_px: int = 30,
    min_pts: int = 6,
    blank_bottom_ticks: bool = True,
) -> list[dict]:
    """Read a black dashed or dotted fit curve at `at` via local linear regression."""
    cL, rT, cR, rB = frame
    mask = rgb_mask(img, (0, 0, 0), frame=frame, inset=6, exclude_fracs=exclude_fracs, black_mode=True)
    rTi, rBi, cLi, cRi = int(round(rT)), int(round(rB)), int(round(cL)), int(round(cR))
    mask[rTi:rTi + 30, :] = False
    if blank_bottom_ticks:
        mask[max(0, rBi - 30):rBi + 1, :] = False
    mask[:, cLi:cLi + 30] = False
    mask[:, max(0, cRi - 30):cRi + 1] = False

    cal_unc_px = max(ya.residual_px, 0.5) + 2.0
    out: list[dict] = []

    for x in at:
        c_nom = xa.to_pixel(x)
        c0 = max(cLi + 6, int(round(c_nom)) - half_win_px)
        c1 = min(cRi - 6, int(round(c_nom)) + half_win_px)
        sub = mask[:, c0:c1 + 1]
        rr, cc_rel = np.nonzero(sub)
        if len(rr) < min_pts:
            out.append({"x": float(x), "y": None, "found": False, "note": "no black curve pixels"})
            continue
        cc = cc_rel.astype(np.float64) + c0 + 0.5
        rr_f = rr.astype(np.float64) + 0.5

        use_branch = branch
        if branch_split_x is not None and x < branch_split_x:
            use_branch = "only"

        order = np.argsort(rr_f)
        rr_s = rr_f[order]
        cc_s = cc[order]
        splits = np.nonzero(np.diff(rr_s) > 35.0)[0] + 1
        r_grps = np.split(rr_s, splits)
        c_grps = np.split(cc_s, splits)
        valid = [(rg, cg) for rg, cg in zip(r_grps, c_grps) if len(rg) >= min_pts]
        if not valid:
            out.append({"x": float(x), "y": None, "found": False, "note": "black fragments too small"})
            continue

        if use_branch == "upper":
            rg, cg = valid[0]
        elif use_branch == "lower":
            rg, cg = valid[-1]
        else:
            rg, cg = max(valid, key=lambda pair: len(pair[0]))

        dx = cg - c_nom
        if np.ptp(dx) >= 4.0:
            slope, intercept = np.polyfit(dx, rg, 1)
            resid = np.abs(rg - (intercept + slope * dx))
            inliers = resid <= max(10.0, 2.5 * float(np.median(resid)) + 3.0)
            if inliers.sum() >= 4 and np.ptp(dx[inliers]) >= 3.0:
                slope, intercept = np.polyfit(dx[inliers], rg[inliers], 1)
            row_y = float(intercept)
            thick_half = max(float(np.percentile(np.abs(rg[inliers] - (intercept + slope * dx[inliers])), 90)), 2.5)
        else:
            row_y = float(np.median(rg))
            thick_half = max(0.5 * float(np.percentile(rg, 90) - np.percentile(rg, 10)), 2.5)

        y_val = ya.to_value(row_y)
        unc_px = thick_half + cal_unc_px
        u1 = ya.to_value(row_y + unc_px)
        u2 = ya.to_value(row_y - unc_px)
        out.append({
            "x": float(x),
            "y": float(y_val),
            "lo": float(min(u1, u2, y_val)),
            "hi": float(max(u1, u2, y_val)),
            "px": [float(c_nom), float(row_y)],
            "px_lo": float(row_y - unc_px),
            "px_hi": float(row_y + unc_px),
            "note": f"black fit ({len(rg)} px)",
        })
    return out


# ---------------------------------------------------------------------------
# Reading Primitive 3: Smooth solid colored curve (e.g. Fig. S1 Eq. S13 curve)
# ---------------------------------------------------------------------------

def read_solid_curve(
    img: np.ndarray,
    xa: D.Axis,
    ya: D.Axis,
    frame: tuple[float, float, float, float],
    rgb: tuple[int, int, int],
    at: list[float],
    *,
    tol: float = 48.0,
    half_win_px: int = 3,
    exclude_fracs: list[tuple[float, float, float, float]] = (),
) -> list[dict]:
    """Read a clean continuous colored curve at `at`."""
    mask = rgb_mask(img, rgb, tol=tol, frame=frame, inset=2, exclude_fracs=exclude_fracs)
    cal_unc_px = max(ya.residual_px, 0.5) + 1.5
    out: list[dict] = []

    for x in at:
        c_nom = xa.to_pixel(x)
        c0 = max(0, int(round(c_nom)) - half_win_px)
        c1 = min(mask.shape[1] - 1, int(round(c_nom)) + half_win_px)
        rr, _ = np.nonzero(mask[:, c0:c1 + 1])
        if len(rr) == 0:
            out.append({"x": float(x), "y": None, "found": False, "note": "curve color not found"})
            continue
        row_y = float(np.median(rr)) + 0.5
        half_h = max(0.5 * float(rr.max() - rr.min() + 1), 2.0)
        y_val = ya.to_value(row_y)
        unc_px = half_h + cal_unc_px
        u1 = ya.to_value(row_y + unc_px)
        u2 = ya.to_value(row_y - unc_px)
        out.append({
            "x": float(x),
            "y": float(y_val),
            "lo": float(min(u1, u2, y_val)),
            "hi": float(max(u1, u2, y_val)),
            "px": [float(c_nom), float(row_y)],
            "px_lo": float(row_y - unc_px),
            "px_hi": float(row_y + unc_px),
            "note": f"solid curve (h={2*half_h:.1f}px)",
        })
    return out


# ---------------------------------------------------------------------------
# Reading Primitive 4: Fig. S1 histogram bars
# ---------------------------------------------------------------------------

def read_figS1_histogram(img: np.ndarray, at: list[float]) -> list[dict]:
    """Calibrate Fig. S1 panel from top x-ticks + left y-ticks and read blue histogram bars."""
    xa, ya, frame = calibrate_matlab_axes(
        img,
        x_all_ticks=[0.0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12],
        y_bounds=(0.0, 25.0),
        y_interior_vals=[5.0, 10.0, 15.0, 20.0],
        x_tick_side="top",
        y_tick_side="left",
    )
    cL, rT, cR, rB = frame
    blue = rgb_mask(img, FIGS1_BAR_BLUE, tol=42.0, frame=frame, inset=2)
    red = rgb_mask(img, MATLAB_ORANGE, tol=55.0, frame=frame, inset=2)
    bin_w_px = abs(xa.to_pixel(0.01) - xa.to_pixel(0.0))
    half_sample = int(0.36 * bin_w_px)
    r0i = int(math.ceil(rT)) + 4
    r1i = int(math.floor(rB)) - 2
    cal_unc_px = max(ya.residual_px, 0.5) + 2.0

    out: list[dict] = []
    for x in at:
        c_nom = xa.to_pixel(x)
        c0 = max(int(math.ceil(cL)) + 3, int(round(c_nom)) - half_sample)
        c1 = min(int(math.floor(cR)) - 3, int(round(c_nom)) + half_sample)
        tops: list[float] = []
        for c in range(c0, c1 + 1):
            ys = np.nonzero(blue[r0i:r1i, c])[0]
            if len(ys) == 0:
                continue
            rt = int(ys[0]) + r0i
            if red[max(r0i, rt - 10):rt, c].any():
                continue
            tops.append(float(rt) - 1.5)
        if not tops:
            any_ys = np.nonzero(blue[r0i:r1i, c0:c1 + 1])[0]
            if len(any_ys) > 0:
                tops = [float(np.percentile(any_ys, 10)) + r0i - 1.5]
            else:
                row_y = float(rB)
                u_hi = ya.to_value(row_y - cal_unc_px)
                out.append({
                    "x": float(x), "y": 0.0, "lo": 0.0, "hi": float(max(0.0, u_hi)),
                    "px": [float(c_nom), row_y], "px_lo": row_y - cal_unc_px, "px_hi": row_y,
                    "note": "zero-height bin",
                })
                continue
        row_y = float(np.percentile(tops, 20))
        y_val = max(0.0, ya.to_value(row_y))
        u1 = max(0.0, ya.to_value(row_y + cal_unc_px))
        u2 = ya.to_value(row_y - cal_unc_px)
        out.append({
            "x": float(x),
            "y": float(y_val),
            "lo": float(min(u1, u2, y_val)),
            "hi": float(max(u1, u2, y_val)),
            "px": [float(c_nom), float(row_y)],
            "px_lo": float(row_y - cal_unc_px),
            "px_hi": float(row_y + cal_unc_px),
            "note": f"hist bar ({len(tops)} clean cols)",
        })
    return out


# ---------------------------------------------------------------------------
# Reading Primitive 5: Fig. 4 bottom-right 3x3 tau_c heatmap via log10 colorbar
# ---------------------------------------------------------------------------

def read_fig4_heatmap_col(img: np.ndarray, kon_rel: float) -> list[dict]:
    """Decode one column (fixed k_on / k_on^(0) in {0.25, 1, 4}) of Fig. 4's 3x3 tau_c heatmap."""
    mx = img.max(axis=2).astype(np.int16)
    mn = img.min(axis=2).astype(np.int16)
    sat = mx - mn
    col_sat = (sat > 30).mean(axis=0) > 0.40
    d = np.diff(np.concatenate(([0], col_sat.view(np.uint8), [0])))
    s = np.nonzero(d == 1)[0]
    e = np.nonzero(d == -1)[0]
    bands = [(int(a), int(b)) for a, b in zip(s, e) if b - a >= 40]
    if len(bands) < 2:
        raise RuntimeError(f"Expected heatmap + colorbar bands, found {bands}")
    hm_c0, hm_c1 = bands[0]
    cb_c0, cb_c1 = bands[1]

    neut_dark = ((mx + mn) < 220) & (sat < 25)
    cb_in_c0, cb_in_c1 = cb_c0 + 8, cb_c1 - 8
    tb_rows = np.nonzero(neut_dark[:, cb_in_c0:cb_in_c1].mean(axis=1) > 0.75)[0]
    if len(tb_rows) < 2:
        raise RuntimeError("Could not find colorbar top/bottom borders")
    r_min, r_max = int(tb_rows[0]), int(tb_rows[-1])
    rT = float(np.mean([r for r in tb_rows if r <= r_min + 6])) + 0.5
    rB = float(np.mean([r for r in tb_rows if r >= r_max - 6])) + 0.5

    tick_strip = neut_dark[int(rT) + 12:int(rB) - 12, max(0, cb_c1 - 18):cb_c1 - 2]
    tick_rows = np.nonzero(tick_strip.sum(axis=1) >= 4)[0] + int(rT) + 12
    grps: list[list[int]] = []
    for r in tick_rows:
        if not grps or r > grps[-1][-1] + 5:
            grps.append([int(r)])
        else:
            grps[-1].append(int(r))
    ticks = [float(np.mean(g)) + 0.5 for g in grps]
    if len(ticks) != 3:
        raise RuntimeError(f"Expected 3 interior colorbar ticks (5, 10, 20), found {len(ticks)}: {ticks}")

    cb_axis = D.Axis("log10", list(zip([rT, ticks[0], ticks[1], ticks[2], rB], [2.0, 5.0, 10.0, 20.0, 60.0])))
    if cb_axis.residual_px > 2.0:
        raise RuntimeError(f"Colorbar log10 fit residual {cb_axis.residual_px:.2f} px")

    cb_r0 = int(math.ceil(rT)) + 2
    cb_r1 = int(math.floor(rB)) - 2
    cb_mid0 = cb_c0 + int(0.25 * (cb_c1 - cb_c0))
    cb_mid1 = cb_c0 + int(0.75 * (cb_c1 - cb_c0))
    cb_rgb = img[cb_r0:cb_r1, cb_mid0:cb_mid1].astype(np.float32).mean(axis=1)

    col_idx = {0.25: 0, 1.0: 1, 4.0: 2}[kon_rel]
    cc_c = hm_c0 + (col_idx + 0.5) * (hm_c1 - hm_c0) / 3.0

    out: list[dict] = []
    for lm, row_i in [(0.2, 2), (0.4, 1), (0.8, 0)]:
        rr_c = rT + (row_i + 0.5) * (rB - rT) / 3.0
        r_ci, c_ci = int(round(rr_c)), int(round(cc_c))
        cell_rgb = img[r_ci - 35:r_ci + 35, c_ci - 35:c_ci + 35].astype(np.float32).mean(axis=(0, 1))
        dists = np.sqrt(((cb_rgb - cell_rgb) ** 2).sum(axis=1))
        dmin = float(dists.min())
        idxs = np.nonzero(dists <= dmin + 2.0)[0]
        row_top = cb_r0 + float(idxs[0]) + 0.5
        row_bot = cb_r0 + float(idxs[-1]) + 0.5
        row_mid = 0.5 * (row_top + row_bot)
        unc_px = max(cb_axis.residual_px, 1.0) + 0.5 * (row_bot - row_top) + 1.5
        val = cb_axis.to_value(row_mid)
        v_a = cb_axis.to_value(max(rT, row_mid - unc_px))
        v_b = cb_axis.to_value(min(rB, row_mid + unc_px))
        out.append({
            "x": float(lm),
            "y": float(val),
            "lo": float(min(v_a, v_b, val)),
            "hi": float(max(v_a, v_b, val)),
            "px": [float(cc_c), float(rr_c)],
            "note": f"heatmap cell (lm={lm}, kon={kon_rel}), cb_row={row_mid:.1f}, dRGB={dmin:.1f}",
        })
    return out


# ---------------------------------------------------------------------------
# Panel-specific convenience wrappers
# ---------------------------------------------------------------------------

def read_figS1_curve(img: np.ndarray, at: list[float]) -> list[dict]:
    xa, ya, frame = calibrate_matlab_axes(
        img,
        x_all_ticks=[0.0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12],
        y_bounds=(0.0, 25.0),
        y_interior_vals=[5.0, 10.0, 15.0, 20.0],
        x_tick_side="top",
        y_tick_side="left",
    )
    return read_solid_curve(img, xa, ya, frame, MATLAB_ORANGE, at, tol=50.0, half_win_px=3)


def read_standard_errorbar(
    img: np.ndarray,
    at: list[float],
    *,
    x_bounds: tuple[float, float],
    y_bounds: tuple[float, float],
    rgb: tuple[int, int, int],
    x_interior: list[float] | None = None,
    y_interior: list[float] | None = None,
    x_scale: str = "linear",
    y_scale: str = "linear",
    exclude_fracs: list[tuple[float, float, float, float]] = (),
    black_mode: bool = False,
    tol: float = 48.0,
    stem_search_px: int = 5,
    expected_y_band: tuple[float, float] | None = None,
    max_ink_gap: int = 35,
    blank_ticks: bool = False,
) -> list[dict]:
    xa, ya, frame = calibrate_matlab_axes(
        img,
        x_bounds=x_bounds,
        y_bounds=y_bounds,
        x_scale=x_scale,
        y_scale=y_scale,
        x_interior_vals=x_interior,
        y_interior_vals=y_interior,
    )
    return read_errorbar_series(
        img, xa, ya, frame, rgb, at,
        tol=tol,
        exclude_fracs=exclude_fracs,
        black_mode=black_mode,
        stem_search_px=stem_search_px,
        expected_y_band=expected_y_band,
        max_ink_gap=max_ink_gap,
        blank_ticks=blank_ticks,
    )


def read_fig2_maxsize(
    img: np.ndarray,
    at: list[float],
    rgb: tuple[int, int, int],
    *,
    tol: float = 48.0,
) -> list[dict]:
    """Fig. 2 bottom-right (log10 y-axis with non-power-of-10 borders ~1.8 to 675)."""
    frame = find_outer_frame(img)
    cL, rT, cR, rB = frame
    xa = D.Axis("linear", [(cL, 0.0), (cR, 60.0)])
    xraw = find_inward_ticks(img, frame, "bottom")
    xmarks = [(cL, 0.0), (cR, 60.0)]
    for v in [10.0, 20.0, 30.0, 40.0, 50.0]:
        pe = xa.to_pixel(v)
        if xraw:
            near = min((p for p, _ in xraw), key=lambda p: abs(p - pe))
            if abs(near - pe) <= 5.0:
                xmarks.append((near, v))
    xa = D.Axis("linear", xmarks)
    ya = calibrate_log_y_from_decades(img, frame, [100.0, 10.0], side="left")
    return read_errorbar_series(img, xa, ya, frame, rgb, at, tol=tol, stem_search_px=5, max_ink_gap=40)


def read_standard_black_fit(
    img: np.ndarray,
    at: list[float],
    *,
    x_bounds: tuple[float, float],
    y_bounds: tuple[float, float],
    x_interior: list[float] | None = None,
    y_interior: list[float] | None = None,
    branch: str = "only",
    branch_split_x: float | None = None,
    exclude_fracs: list[tuple[float, float, float, float]] = (),
    half_win_px: int = 32,
) -> list[dict]:
    xa, ya, frame = calibrate_matlab_axes(
        img,
        x_bounds=x_bounds,
        y_bounds=y_bounds,
        x_interior_vals=x_interior,
        y_interior_vals=y_interior,
    )
    return read_black_fit_curve(
        img, xa, ya, frame, at,
        branch=branch,
        branch_split_x=branch_split_x,
        exclude_fracs=exclude_fracs,
        half_win_px=half_win_px,
    )


def read_figS3_hist(
    img: np.ndarray,
    at: list[float],
    *,
    x_bounds: tuple[float, float],
    y_bounds: tuple[float, float],
    rgb: tuple[int, int, int],
    x_scale: str = "linear",
    x_interior: list[float] | None = None,
    y_interior: list[float] | None = None,
    exclude_fracs: list[tuple[float, float, float, float]] = (),
    tol: float = 48.0,
) -> list[dict]:
    """Read discrete hollow-circle markers + vertical error bars in Fig. S3.

    Uses a 4-quadrant diagonal ring template (radius 20-25 px, |dx|>=10, |dy|>=10)
    to pinpoint the exact circle marker center (immune to horizontal error-bar
    caps, vertical stems, and bottom clipping at y=0), then walks up/down from
    the marker center along the vertical stem to find bar_hi and bar_lo.
    """
    xa, ya, frame = calibrate_matlab_axes(
        img,
        x_bounds=x_bounds,
        y_bounds=y_bounds,
        x_scale=x_scale,
        y_scale="linear",
        x_interior_vals=x_interior,
        y_interior_vals=y_interior,
        x_tick_side="top",
        y_tick_side="left",
    )
    cL, rT, cR, rB = frame
    H, W = img.shape[:2]
    ext_frame = (max(0.0, cL - 40.0), max(0.0, rT - 40.0), min(float(W), cR + 40.0), min(float(H), rB + 40.0))
    mask = rgb_mask(img, rgb, tol=tol, frame=ext_frame, inset=0)
    for fx0, fy0, fx1, fy1 in exclude_fracs:
        ec0 = int(round(cL + fx0 * (cR - cL)))
        ec1 = int(round(cL + fx1 * (cR - cL)))
        er0 = int(round(rT + fy0 * (rB - rT)))
        er1 = int(round(rT + fy1 * (rB - rT)))
        mask[max(0, er0):min(H, er1), max(0, ec0):min(W, ec1)] = False

    other_ink = (img.min(axis=2) < 225) & (~mask)
    cal_unc_px = max(ya.residual_px, 0.5) + 2.0

    # Precompute 4-quadrant diagonal ring template (excludes stem |dx|<10 and cap |dy|<10)
    ang = np.linspace(0.0, 2.0 * np.pi, 32, endpoint=False)
    pts = []
    for R in (20.0, 25.0):
        dx = np.round(R * np.cos(ang)).astype(int)
        dy = np.round(R * np.sin(ang)).astype(int)
        keep = (np.abs(dx) >= 10) & (np.abs(dy) >= 10)
        for a, b in zip(dx[keep], dy[keep]):
            pts.append((int(a), int(b)))
    dxs = np.array([p[0] for p in pts], dtype=int)
    dys = np.array([p[1] for p in pts], dtype=int)
    rs_scan = np.arange(int(round(rT)), int(round(rB)) + 1, dtype=int)

    out: list[dict] = []
    for x in at:
        c_nom = xa.to_pixel(x)
        c_int = int(round(c_nom))

        # 1. Detect circle marker center via 4-quadrant ring score
        rr = rs_scan[:, None] + dys[None, :]
        cc = c_int + dxs[None, :]
        ok = (rr >= 0) & (rr < H) & (cc >= 0) & (cc < W)
        rr_c = np.clip(rr, 0, H - 1)
        cc_c = np.clip(cc, 0, W - 1)
        hits = (mask[rr_c, cc_c] & ok).sum(axis=1)
        denom = np.maximum(1, ok.sum(axis=1))
        scores = hits / denom
        mx_score = float(scores.max())
        r_pk = int(rs_scan[int(np.argmax(scores))])
        near = rs_scan[(scores >= mx_score - 0.04) & (np.abs(rs_scan - r_pk) <= 14)]
        row_marker = float(near.mean()) + 0.5

        if mx_score < 0.30:
            # Fallback for very small/partially occluded markers: check vertical stem
            out.append({"x": float(x), "y": None, "found": False, "note": f"circle occluded (score={mx_score:.2f})"})
            continue

        # Clamp marker very close to rB (y=0) cleanly to 0
        if abs(row_marker - rB) <= 12.0:
            y_val = 0.0
            row_y = float(rB)
        else:
            y_val = max(0.0, ya.to_value(row_marker))
            row_y = ya.to_pixel(y_val)

        # 2. Trace errorbar stem upward & downward from row_y in [c_int - 5 .. c_int + 5]
        # Walk up/down bridging another curve's ink (up to 65 px) and tiny 2-px gaps
        c_lo_s = max(0, c_int - 5)
        c_hi_s = min(W - 1, c_int + 5)
        col_m = mask[:, c_lo_s:c_hi_s + 1].any(axis=1)
        col_o = other_ink[:, c_lo_s:c_hi_s + 1].any(axis=1)
        r_cur = int(round(row_y))

        # Walk up from r_cur
        r_top = r_cur
        r_scan = r_cur - 1
        while r_scan >= int(round(rT)) - 4:
            if col_m[r_scan]:
                r_top = r_scan
                r_scan -= 1
            else:
                # Look ahead up to 65 px if covered by other_ink (or 3 px white)
                found_next = None
                for k in range(1, 66):
                    rn = r_scan - k
                    if rn < int(round(rT)) - 4:
                        break
                    if col_m[rn]:
                        if k <= 3 or col_o[rn + 1:r_scan + 1].mean() >= 0.70:
                            found_next = rn
                        break
                if found_next is not None:
                    r_top = found_next
                    r_scan = found_next - 1
                else:
                    break

        # Walk down from r_cur
        r_bot = r_cur
        r_scan = r_cur + 1
        while r_scan <= int(round(rB)) + 4:
            if col_m[r_scan]:
                r_bot = r_scan
                r_scan += 1
            else:
                found_next = None
                for k in range(1, 66):
                    rn = r_scan + k
                    if rn > int(round(rB)) + 4:
                        break
                    if col_m[rn]:
                        if k <= 3 or col_o[r_scan:rn].mean() >= 0.70:
                            found_next = rn
                        break
                if found_next is not None:
                    r_bot = found_next
                    r_scan = found_next + 1
                else:
                    break

        # Note: marker outer circle extends ~31 px above/below center; if no errorbar protrudes
        # beyond the circle radius (~34 px), the errorbar is smaller than the marker (or zero).
        # Otherwise r_top / r_bot mark the exact cap outer edges!
        clipped_bot = r_bot >= int(round(rB)) - 3
        clipped_top = r_top <= int(round(rT)) + 3
        if y_val == 0.0:
            r_hi_px = float(r_top) + 3.0 if (rB - r_top) > 36 else float(rB)
            r_lo_px = float(rB)
        else:
            r_hi_px = float(r_top) + (0.5 if clipped_top else 3.5) if (row_y - r_top) > 34 else row_y
            r_lo_px = float(min(rB, r_bot)) + (0.5 if clipped_bot else -2.5) if (r_bot - row_y) > 34 else row_y
            # If one side of the symmetric errorbar was clipped/occluded by another marker,
            # verify using symmetry where valid (only if one bar clearly protruded and the other didn't
            # ONLY because it hit rB=0 -> keep clipped at 0!).
            if r_hi_px > row_y:
                r_hi_px = row_y
            if r_lo_px < row_y:
                r_lo_px = row_y

        bar_hi = max(y_val, ya.to_value(r_hi_px))
        bar_lo = max(0.0, min(y_val, ya.to_value(r_lo_px)))
        u1 = max(0.0, ya.to_value(row_y + cal_unc_px))
        u2 = ya.to_value(row_y - cal_unc_px)
        bh_out = max(bar_hi, ya.to_value(r_hi_px - cal_unc_px))
        bl_out = max(0.0, min(bar_lo, ya.to_value(r_lo_px + cal_unc_px)))

        out.append({
            "x": float(x),
            "y": float(y_val),
            "lo": float(min(u1, u2, y_val)),
            "hi": float(max(u1, u2, y_val)),
            "bar_lo": float(bar_lo),
            "bar_hi": float(bar_hi),
            "bar_lo_outer": float(bl_out),
            "bar_hi_outer": float(bh_out),
            "px": [float(c_nom), float(row_y)],
            "px_lo": float(r_hi_px),
            "px_hi": float(r_lo_px),
            "note": f"circle score={mx_score:.2f}, bar=[{bar_lo:.2f},{bar_hi:.2f}]",
        })
    return out


