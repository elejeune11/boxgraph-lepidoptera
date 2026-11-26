"""Calibrated image reader for S12 Figure panel (b) (p. 37 / SI-18): pairwise filament-end distance
probability density curves Dis--, Dis+-, and Dis++ (nm^-1 vs nm, 0..8000 nm, y: 0..0.01) for two
unipolar non-treadmilling conditions at alpha:A = 0.05: M:A = 0.0225 (left column, bundle-like) and
M:A = 0.675 (right column, aster-like) (1200 dpi).

Detects the three vertically stacked 4-spine MATLAB subpanels in each column crop, calibrates each
box independently on its 9 x-ticks (0..8000 nm) and 3 y-ticks (0, 0.005, 0.01), and extracts the
solid black mean probability density curve embedded inside the red scatter envelope.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _common as C


def find_s12b_subpanels(image: np.ndarray) -> list[dict]:
    """Locate the 3 vertically stacked 4-spine boxes (Dis--, Dis+-, Dis++) in one S12(b) column."""
    r = image[:, :, 0]
    g = image[:, :, 1]
    b = image[:, :, 2]
    gray3 = r.astype(np.uint16) + g + b
    cmax = np.maximum(np.maximum(r, g), b)
    cmin = np.minimum(np.minimum(r, g), b)
    # Brightness-only mask for spines (red cloud at y=0 sits on top of the black bottom spine)
    dark_sp = gray3 < 540
    # Neutral dark mask for ticks
    dark_tk = (gray3 < 540) & ((cmax - cmin) < 60)
    H, W = dark_sp.shape

    h_sp = C.find_spines_1d(dark_sp, min_len=0.55 * W, gap=6, merge_dist=5)
    # Exclude any stray border line within 15 px of the crop top/bottom edges
    h_sp = [s for s in h_sp if 15 < s["lo"] < H - 15]
    if len(h_sp) < 6:
        raise ValueError(f"S12(b): expected 6 horizontal spines for 3 stacked boxes, got {len(h_sp)}")
    # Keep the 6 spines with matching horizontal box span (>= 85% of median of top 6)
    by_len = sorted(h_sp, key=lambda s: -s["length"])[:6]
    h_sp = sorted(by_len, key=lambda s: s["lo"])

    panels = []
    for k in range(3):
        top_sp, bot_sp = h_sp[2 * k], h_sp[2 * k + 1]
        sub_d = dark_tk[top_sp["lo"]:bot_sp["hi"] + 1, :]
        sub_H = sub_d.shape[0]
        v_sp = C.find_spines_1d(sub_d.T, min_len=0.55 * sub_H, gap=5, merge_dist=5)
        if len(v_sp) < 2:
            raise ValueError(f"S12(b) row {k}: expected left & right spines, got {len(v_sp)}")
        left_sp = min(v_sp, key=lambda s: abs(s["center"] - top_sp["start"]))
        right_sp = min(v_sp, key=lambda s: abs(s["center"] - top_sp["end"]))
        c0, c1 = left_sp["hi"] + 1, right_sp["lo"]
        r0, r1 = top_sp["hi"] + 1, bot_sp["lo"]

        depth = 28
        x_top_strip = dark_tk[top_sp["hi"] + 1:min(H, top_sp["hi"] + 1 + depth), :]
        xt_top = C.find_ticks_in_strip(x_top_strip, c0 + 8, c1 - 8, min_len=4, max_len=24, max_width=16)
        x_bot_strip = dark_tk[max(0, bot_sp["lo"] - depth):bot_sp["lo"], :][::-1, :]
        xt_bot = C.find_ticks_in_strip(x_bot_strip, c0 + 8, c1 - 8, min_len=4, max_len=24, max_width=16)

        y_left_strip = dark_tk[:, left_sp["hi"] + 1:min(W, left_sp["hi"] + 1 + depth)].T
        yt_left = C.find_ticks_in_strip(y_left_strip, r0 + 8, r1 - 8, min_len=4, max_len=24, max_width=16)
        y_right_strip = dark_tk[:, max(0, right_sp["lo"] - depth):right_sp["lo"]][:, ::-1].T
        yt_right = C.find_ticks_in_strip(y_right_strip, r0 + 8, r1 - 8, min_len=4, max_len=24, max_width=16)

        # Match detected inner x-ticks against uniform 0..8000 nm grid spanning [left_sp, right_sp]
        x_vals_all = [0.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.0]
        exp_x = np.linspace(left_sp["center"], right_sp["center"], 9)
        x_px, x_v = [left_sp["center"]], [0.0]
        for ex, xv in zip(exp_x[1:-1], x_vals_all[1:-1]):
            near = [t for t in (xt_top + xt_bot) if abs(t - ex) <= 14.0]
            if near:
                x_px.append(min(near, key=lambda t: abs(t - ex)))
                x_v.append(xv)
        x_px.append(right_sp["center"])
        x_v.append(8000.0)

        yt = C.merge_tick_set(top_sp["center"], yt_right + yt_left, bot_sp["center"], 3, tol_px=18.0)
        xa = C.fit_axis(x_px, x_v, max_residual_px=3.5)
        ya = C.fit_axis(yt, [0.01, 0.005, 0.0], max_residual_px=3.5)
        panels.append({
            "box": (int(c0), int(r0), int(c1), int(r1)),
            "bottom": bot_sp["center"],
            "xa": xa,
            "ya": ya,
        })
    return panels


def read_s12b_curve(image: np.ndarray, row_idx: int, label: str, at: list[float]) -> list[dict]:
    """Read the solid black mean density curve from S12(b) subpanel row_idx (0=Dis--, 1=Dis+-, 2=Dis++)."""
    panels = find_s12b_subpanels(image)
    p = panels[row_idx]
    c0, r0, c1, r1 = p["box"]
    xa, ya = p["xa"], p["ya"]
    bot_center = p["bottom"]

    r = image[:, :, 0].astype(np.float64)
    g = image[:, :, 1].astype(np.float64)
    b = image[:, :, 2].astype(np.float64)
    red_cloud = (r > 165) & ((r - g) > 85) & ((r - b) > 85)

    out = []
    for xv in at:
        pc = xa.to_pixel(xv)
        col = int(round(pc))
        ca = max(c0 + 2, col - 4)
        cb = min(c1 - 2, col + 5)
        y_top_lim = r0 + 18
        y_bot_lim = r1 - 5
        sub_red = red_cloud[y_top_lim:y_bot_lim, ca:cb]

        if not sub_red.any():
            row_val = bot_center
            half_px = 2.0
            tag = "baseline y=0"
        else:
            prof_r = r[y_top_lim:y_bot_lim, ca:cb].mean(axis=1)
            prof_g = g[y_top_lim:y_bot_lim, ca:cb].mean(axis=1)
            prof_b = b[y_top_lim:y_bot_lim, ca:cb].mean(axis=1)
            # Continuous black line (possibly red-blended) has low G, low B, and strictly lower R than bright red dots
            cand = (prof_g < 55.0) & (prof_b < 55.0) & (prof_r < 195.0)
            if not cand.any():
                # Fallback for very narrow steep flank where ±4 px window spans white background: use ±2 px window
                ca2, cb2 = max(c0 + 2, col - 2), min(c1 - 2, col + 3)
                prof_r = r[y_top_lim:y_bot_lim, ca2:cb2].mean(axis=1)
                prof_g = g[y_top_lim:y_bot_lim, ca2:cb2].mean(axis=1)
                prof_b = b[y_top_lim:y_bot_lim, ca2:cb2].mean(axis=1)
                cand = (prof_g < 85.0) & (prof_b < 85.0) & (prof_r < 145.0)
            if not cand.any():
                row_val = bot_center
                half_px = 2.0
                tag = "baseline y=0"
            else:
                score = np.where(cand, prof_r + 0.5 * prof_g, 999.0)
                k = int(np.argmin(score))
                best_row = y_top_lim + k
                # Refine centroid over nearby rows within +-5 px whose score is within +22 of minimum
                lo_k = max(0, k - 5)
                hi_k = min(len(score), k + 6)
                near = np.flatnonzero(score[lo_k:hi_k] <= score[k] + 22.0) + (y_top_lim + lo_k)
                row_val = 0.5 * (float(near.min()) + float(near.max())) + 0.5
                half_px = max(2.0, 0.5 * (float(near.max()) - float(near.min()) + 1.0))
                tag = "black mean curve"

        slope = 0.0
        if col - 5 >= c0 + 2 and col + 5 < c1 - 2 and tag != "baseline y=0":
            # Quick slope from minimum-R ridge at col-5 and col+5 near row_val
            ra = max(y_top_lim, int(round(row_val)) - 20)
            rb = min(y_bot_lim, int(round(row_val)) + 21)
            if rb > ra + 4:
                pl = r[ra:rb, col - 6:col - 3].mean(axis=1) + 0.5 * g[ra:rb, col - 6:col - 3].mean(axis=1)
                pr = r[ra:rb, col + 4:col + 7].mean(axis=1) + 0.5 * g[ra:rb, col + 4:col + 7].mean(axis=1)
                slope = (float(np.argmin(pr)) - float(np.argmin(pl))) / 10.0

        note = f"S12(b) {label} at d={xv:g} nm ({tag}; rx={xa.residual_px:.2f}px, ry={ya.residual_px:.2f}px)"
        out.append(
            C.make_reading(
                xv, row_val, ya, pc, half_px=half_px, x_axis=xa, slope_px=slope,
                note=note, clamp_y=(0.0, 0.01),
            )
        )
    return out
