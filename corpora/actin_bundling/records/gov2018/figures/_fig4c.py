"""Shared calibrated reader for Fig. 4(c) (p. 8): omega'(q) and omega''(q) dispersion relations."""

from __future__ import annotations

import numpy as np
from _calib import (
    clusters,
    find_inward_ticks,
    find_neutral_frame,
    match_ticks_regular,
    runs_1d,
    sample_dashed_curve_y,
    sample_solid_curve_y,
    verify_axis,
)


def _calibrate_fig4c(image: np.ndarray):
    dark, gray, h_lines, v_lines = find_neutral_frame(
        image, row_frac=0.50, col_frac=0.30, gray_thresh=125.0, sat_thresh=35
    )
    if len(v_lines) < 2 or len(h_lines) < 5:
        raise RuntimeError(f"Fig. 4(c): unexpected frame lines h={h_lines}, v={v_lines}")

    v_left, v_right = v_lines[0], v_lines[-1]
    # h_lines are: top_frame_top, top_zero_line, top_frame_bot, bot_frame_top, bot_frame_bot
    h_top0, h_zero, h_top1, h_bot0, h_bot1 = h_lines[0], h_lines[1], h_lines[2], h_lines[3], h_lines[4]

    c0, c1 = int(round(v_left[0])), int(round(v_right[0]))

    # --- Top panel x-axis: q = 0, 2, 4, ..., 20 ---
    _, xt_top_b, yt_top_l, _ = find_inward_ticks(
        dark, int(round(h_top0[0])), int(round(h_top1[0])), c0, c1, band=(6, 18), margin=15, max_width=15, frac=0.55
    )
    x_cands_top = [v_left[0]] + xt_top_b + [v_right[0]]
    x_ticks_top = match_ticks_regular(x_cands_top, 11, tol=6.0)
    x_top = verify_axis("linear", list(zip(x_ticks_top, [2.0 * k for k in range(11)])), max_residual_px=4.0)

    # --- Top panel y-axis: omega' = 0.02, 0, -0.02, -0.04, -0.06, -0.08, -0.10, -0.12 ---
    # Left border ticks are unobstructed at all 8 values
    y_ticks_top = match_ticks_regular(yt_top_l, 8, tol=8.0)
    y_top_vals = [0.02, 0.0, -0.02, -0.04, -0.06, -0.08, -0.10, -0.12]
    y_top = verify_axis("linear", list(zip(y_ticks_top, y_top_vals)), max_residual_px=5.0)

    # --- Bottom panel x-axis: q = 0, 2, 4, ..., 20 ---
    _, xt_bot_b, yt_bot_l, yt_bot_r = find_inward_ticks(
        dark, int(round(h_bot0[0])), int(round(h_bot1[0])), c0, c1, band=(6, 18), margin=15, max_width=15, frac=0.55
    )
    x_cands_bot = [v_left[0]] + xt_bot_b + [v_right[0]]
    x_ticks_bot = match_ticks_regular(x_cands_bot, 11, tol=6.0)
    x_bot = verify_axis("linear", list(zip(x_ticks_bot, [2.0 * k for k in range(11)])), max_residual_px=4.0)

    # --- Bottom panel left y-axis: 0.15 (top border), 0.10, 0.05 (inward ticks), 0.00 (extrapolated / red baseline) ---
    # yt_bot_l gives the two interior ticks [0.10, 0.05]; top border h_bot0[0] is 0.15
    yt_l_sorted = sorted(yt_bot_l)
    if len(yt_l_sorted) < 2:
        raise RuntimeError(f"Fig. 4(c) bottom left y-ticks missing: {yt_l_sorted}")
    # Also detect the exact baseline row of red/green at q = 1..8 (left zero line)
    r_ch = image[:, :, 0].astype(np.int16)
    g_ch = image[:, :, 1].astype(np.int16)
    b_ch = image[:, :, 2].astype(np.int16)
    red_base = (r_ch > 180) & (r_ch - g_ch > 60) & (r_ch - b_ch > 60)
    rb_rows, _ = np.nonzero(red_base[int(h_bot1[0]) - 80 : int(h_bot1[0]) - 15, c0 + 100 : c0 + 800])
    y0_left_px = float(np.median(rb_rows)) + (int(h_bot1[0]) - 80) + 0.5
    # Verify spacing: top border (0.15) -> yt_l[0] (0.10) -> yt_l[1] (0.05) -> y0_left_px (0.00)
    y_bot_left = verify_axis(
        "linear",
        [(h_bot0[0], 0.15), (yt_l_sorted[0], 0.10), (yt_l_sorted[1], 0.05), (y0_left_px, 0.0)],
        max_residual_px=3.0,
    )

    # --- Bottom panel right y-axis: 5 (top border), 4, 3, 2, 1 (inward ticks) ---
    yt_r_4 = match_ticks_regular(yt_bot_r, 4, tol=4.0)
    y_bot_right = verify_axis(
        "linear",
        [(h_bot0[0], 5.0), (yt_r_4[0], 4.0), (yt_r_4[1], 3.0), (yt_r_4[2], 2.0), (yt_r_4[3], 1.0)],
        max_residual_px=3.0,
    )

    boxes = {
        "top": (int(h_top0[2]) + 2, int(h_top1[1]) - 2, int(v_left[2]) + 2, int(v_right[1]) - 2),
        "bot": (int(h_bot0[2]) + 2, int(h_bot1[1]) - 2, int(v_left[2]) + 2, int(v_right[1]) - 2),
    }
    return x_top, y_top, x_bot, y_bot_left, y_bot_right, boxes


def _color_masks(image: np.ndarray) -> dict[str, np.ndarray]:
    r = image[:, :, 0].astype(np.int16)
    g = image[:, :, 1].astype(np.int16)
    b = image[:, :, 2].astype(np.int16)
    return {
        "cyan": (g > 155) & (b > 150) & (g - r > 22) & (b - r > 15) & (np.abs(g - b) < 35),
        "green": (g > 110) & (g - r > 20) & (g - b > 35),
        "red": (r > 180) & (r - g > 60) & (r - b > 60),
        "darkblue": (b > 130) & (b - r > 35) & (b - g > 20),
        "purple": (r > 130) & (b > 130) & (r - g > 28) & (b - g > 28) & (np.abs(r - b) < 30),
    }


def _trace_upper_cyan(mask: np.ndarray, rmin: int, rmax: int, cmin: int, cmax: int) -> np.ndarray:
    """Return column-wise row center of the solid (upper) cyan curve across [cmin..cmax]."""
    W = mask.shape[1]
    upper = np.full(W, np.nan)
    for c in range(cmin, cmax + 1):
        rs = runs_1d(mask[rmin : rmax + 1, c], gap=2)
        if rs:
            a, b = rs[0]  # topmost cyan run in this column
            upper[c] = 0.5 * (a + b) + rmin + 0.5
    # Fill tiny gaps by linear interpolation
    good = np.isfinite(upper)
    if good.sum() >= 2:
        idx = np.arange(W)
        upper = np.interp(idx, idx[good], upper[good])
    return upper


def read_fig4c(image: np.ndarray, at: list[float], curve: str) -> list[dict]:
    x_top, y_top, x_bot, y_bot_left, y_bot_right, boxes = _calibrate_fig4c(image)
    masks = _color_masks(image)

    if curve.startswith("omega_real_"):
        rmin, rmax, cmin, cmax = boxes["top"]
        tag = curve[len("omega_real_") :]
        if tag == "cyan":
            # Solid light-blue curve is the topmost cyan curve for all q
            return [
                sample_solid_curve_y(
                    masks["cyan"], x_top, y_top, x, rmin, rmax,
                    pick_fn=lambda rs, c: rs[0],
                    note_prefix="fig4c top cyan solid; ",
                )
                for x in at
            ]
        if tag == "nonlocal":
            # Dashed light-blue curve lies strictly below the solid cyan curve for q >= 1.2
            upper = _trace_upper_cyan(masks["cyan"], rmin, rmax, cmin, cmax)
            dash_mask = masks["cyan"].copy()
            H, W = dash_mask.shape
            rows_grid = np.arange(H)[:, None]
            # Exclude pixels at or above (upper + 18 px)
            for c in range(cmin, cmax + 1):
                if np.isfinite(upper[c]):
                    cut = int(round(upper[c] + 18.0))
                    dash_mask[: min(H, max(0, cut)), c] = False
            return [
                sample_dashed_curve_y(
                    dash_mask, x_top, y_top, x, rmin, rmax, max_search_dx=95,
                    note_prefix="fig4c top cyan dashed; ",
                )
                for x in at
            ]
        if tag == "green":
            return [
                sample_solid_curve_y(
                    masks["green"], x_top, y_top, x, rmin, rmax,
                    note_prefix="fig4c top green; ",
                )
                for x in at
            ]
        if tag == "red":
            return [
                sample_solid_curve_y(
                    masks["red"], x_top, y_top, x, rmin, rmax,
                    search_dx=8,  # bridge small cyan-dash overlaps on q >= 15
                    note_prefix="fig4c top red; ",
                )
                for x in at
            ]
        if tag == "blue":
            return [
                sample_solid_curve_y(
                    masks["darkblue"], x_top, y_top, x, rmin, rmax,
                    note_prefix="fig4c top dark blue (A*=0); ",
                )
                for x in at
            ]

    if curve.startswith("omega_imag_"):
        rmin, rmax, cmin, cmax = boxes["bot"]
        tag = curve[len("omega_imag_") :]
        if tag == "purple":
            return [
                sample_solid_curve_y(
                    masks["purple"], x_bot, y_bot_right, x, rmin, rmax,
                    note_prefix="fig4c bot purple (right axis 0..5 Hz); ",
                )
                for x in at
            ]
        if tag == "cyan":
            # Fit purple straight line so we know where solid cyan sits (just below purple, never down on dashed cyan)
            pr, pc = np.nonzero(masks["purple"][rmin : rmax + 1, cmin : cmax + 1])
            p_m, p_b = np.polyfit(pc + cmin + 0.5, pr + rmin + 0.5, 1)

            def pick_solid_cyan(rs, c):
                # Solid cyan lies just below purple (within [p(c) - 10, p(c) + 120]), well above dashed cyan
                p_row = p_m * (c + 0.5) + p_b
                cands = [seg for seg in rs if (p_row - 10.0) <= 0.5 * (seg[0] + seg[1]) <= (p_row + 120.0)]
                return cands[0] if cands else None

            return [
                sample_solid_curve_y(
                    masks["cyan"], x_bot, y_bot_right, x, rmin, rmax,
                    pick_fn=pick_solid_cyan,
                    search_dx=8,
                    note_prefix="fig4c bot cyan solid (right axis 0..5 Hz); ",
                )
                for x in at
            ]
        if tag == "nonlocal":
            # Dashed light-blue curve lies strictly below the purple line + solid cyan band
            pr, pc = np.nonzero(masks["purple"][rmin : rmax + 1, cmin : cmax + 1])
            p_m, p_b = np.polyfit(pc + cmin + 0.5, pr + rmin + 0.5, 1)
            upper = _trace_upper_cyan(masks["cyan"], rmin, rmax, cmin, cmax)
            dash_mask = masks["cyan"].copy()
            H = dash_mask.shape[0]
            for c in range(cmin, cmax + 1):
                p_row = p_m * (c + 0.5) + p_b
                u_row = upper[c] if (np.isfinite(upper[c]) and abs(upper[c] - p_row) < 120.0) else p_row
                cut = int(round(max(u_row, p_row) + 18.0))
                dash_mask[: min(H, max(0, cut)), c] = False
            return [
                sample_dashed_curve_y(
                    dash_mask, x_bot, y_bot_right, x, rmin, rmax, max_search_dx=95,
                    note_prefix="fig4c bot cyan dashed (right axis 0..5 Hz); ",
                )
                for x in at
            ]
        if tag == "red":
            return [
                sample_solid_curve_y(
                    masks["red"], x_bot, y_bot_left, x, rmin, rmax,
                    note_prefix="fig4c bot red (left axis 0..0.15 Hz); ",
                )
                for x in at
            ]
        if tag == "green":
            # Green curve lies on left axis 0..0.15 Hz with small complex arch near q ~ 18.65
            return [
                sample_solid_curve_y(
                    masks["green"], x_bot, y_bot_left, x, rmin, rmax,
                    pick_fn=lambda rs, c: rs[0],  # topmost green pixel (arch top over zero line)
                    note_prefix="fig4c bot green (left axis 0..0.15 Hz); ",
                )
                for x in at
            ]

    raise ValueError(f"Unknown Fig. 4(c) curve: {curve}")
