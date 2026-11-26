"""Calibrated image reader for S11 Figure (p. 36 / SI-17): treadmilling rate (monomers/(fil. s))
vs treadmilling scale factor chi across 6 panels — Unipolar (a-c) and Apolar (d-f) at
alpha:A = 0.01, 0.1, 0.4 — and 4 M:A ratios (0.0225 red, 0.09 blue, 0.225 green, 0.675 black) (900 dpi).

Calibrates each 4-spine MATLAB panel independently using top/bottom x-ticks (0, 5, 10) and
left/right y-ticks (0..7), isolates each series via orthogonal RGB channel masks, recovers exact
data vertices by fitting the straight connecting line segments between consecutive chi nodes,
and extracts symmetric +/- SD error-bar caps in each series's own color mask.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _common as C

CHI_BLACK = [0.3, 0.6, 1.0, 3.0, 6.0, 10.0]
CHI_COLOR = [1.0, 3.0, 6.0, 10.0]


def calibrate_s11(image: np.ndarray):
    """Find the 4-spine plot box and calibrate x (0, 5, 10) and y (0..7) for one S11 panel."""
    info = C.find_matlab_box(
        image, dark_thresh=160.0, min_h_frac=0.50, min_v_frac=0.50,
        has_top=True, has_right=True, min_tick_len=5, max_tick_len=45,
        max_tick_w=24, margin_px=10, neutral_only=False,
    )
    sp = info["spines"]
    info_n = C.find_matlab_box(
        image, dark_thresh=160.0, min_h_frac=0.50, min_v_frac=0.50,
        has_top=True, has_right=False, min_tick_len=5, max_tick_len=45,
        max_tick_w=24, margin_px=10, neutral_only=True,
    )
    xt = C.merge_tick_set(sp["left"], info_n["x_ticks_top"] + info_n["x_ticks_bot"], sp["right"], 3, tol_px=20.0)
    yt = C.merge_tick_set(sp["top"], info_n["y_ticks_left"], sp["bottom"], 8, tol_px=18.0)
    xa = C.fit_axis(xt, [0.0, 5.0, 10.0], max_residual_px=2.5)
    ya = C.fit_axis(yt, [7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0], max_residual_px=2.5)
    return info, xa, ya


def _raw_color_mask(image: np.ndarray, color: str, box: tuple[int, int, int, int],
                    has_legend: bool = False) -> np.ndarray:
    r = image[:, :, 0].astype(np.int16)
    g = image[:, :, 1].astype(np.int16)
    b = image[:, :, 2].astype(np.int16)
    if color == "red":
        m = (r > 160) & (g < 105) & (b < 105) & ((r - g) > 80) & ((r - b) > 80)
    elif color == "blue":
        m = (b > 160) & (r < 115) & (g < 125) & ((b - r) > 70) & ((b - g) > 60)
    elif color == "green":
        m = (g > 155) & (b < 100) & ((g - b) > 95) & ((g - r) > 40)
    elif color == "black":
        cmax = np.maximum(np.maximum(r, g), b)
        cmin = np.minimum(np.minimum(r, g), b)
        m = (r + g + b < 320) & ((cmax - cmin) < 42)
    else:
        raise ValueError(color)

    if has_legend:
        c0, r0, c1, r1 = box
        # Legend box sits in upper-right interior, stopping cleanly ~35 px LEFT of the right spine
        m[r0 + 2:r0 + int(0.36 * (r1 - r0)), c0 + int(0.54 * (c1 - c0)):c1 - 30] = False
    return m


def _line_mask(mask_raw: np.ndarray, color: str, box: tuple[int, int, int, int]) -> np.ndarray:
    c0, r0, c1, r1 = box
    H, W = mask_raw.shape
    out = np.zeros_like(mask_raw)
    if color == "black":
        out[r0 + 8:r1 - 6, c0 + 6:c1 - 6] = mask_raw[r0 + 8:r1 - 6, c0 + 6:c1 - 6]
        # Mask panel letter "a)".."f)" and inward ticks on top/bottom/left spines
        out[r0:r0 + int(0.15 * (r1 - r0)), c0:c0 + int(0.22 * (c1 - c0))] = False
        out[r0:r0 + 24, :] = False
        out[r1 - 24:r1 + 1, c0 + int(0.45 * (c1 - c0)):c0 + int(0.55 * (c1 - c0))] = False
        out[:, c0:c0 + 22] = False
    else:
        out[r0 + 6:r1 - 3, c0 + 3:min(W, c1 + 22)] = mask_raw[r0 + 6:r1 - 3, c0 + 3:min(W, c1 + 22)]
    return out


def _fit_line_between(mask: np.ndarray, col_L: float, col_R: float,
                      pad_L: int = 24, pad_R: int = 24) -> tuple[float, float, float] | None:
    """Fit the single straight connecting line of this color in [col_L + pad_L, col_R - pad_R]."""
    c_a = int(round(col_L)) + pad_L
    c_b = int(round(col_R)) - pad_R
    if c_b - c_a < 10:
        c_a = int(round(col_L)) + max(6, int(0.2 * (col_R - col_L)))
        c_b = int(round(col_R)) - max(6, int(0.2 * (col_R - col_L)))
    if c_b <= c_a:
        return None
    cols = np.arange(c_a, c_b + 1, 2)
    xs, ys = [], []
    for c in cols:
        rs = np.flatnonzero(mask[:, c])
        if rs.size == 0:
            continue
        sp = np.flatnonzero(np.diff(rs) > 3)
        st = np.concatenate(([0], sp + 1))
        en = np.concatenate((sp, [rs.size - 1]))
        runs = [(0.5 * (rs[a] + rs[b]) + 0.5, int(rs[b] - rs[a] + 1)) for a, b in zip(st, en) if 4 <= (rs[b] - rs[a] + 1) <= 32]
        if len(runs) == 1:
            xs.append(float(c) + 0.5)
            ys.append(runs[0][0])
        elif len(runs) > 1:
            # Prefer thickest valid line run (~8-12 px), rejecting thin anti-aliasing overlap speckles
            runs.sort(key=lambda t: -t[1])
            xs.append(float(c) + 0.5)
            ys.append(runs[0][0])
    if len(xs) < 4:
        return None
    xa_arr = np.array(xs, dtype=float)
    ya_arr = np.array(ys, dtype=float)
    m, b = np.polyfit(xa_arr, ya_arr, 1)
    res = np.abs(ya_arr - (m * xa_arr + b))
    inliers = res <= max(3.0, float(np.percentile(res, 65)) * 2.0)
    if inliers.sum() >= 3:
        m, b = np.polyfit(xa_arr[inliers], ya_arr[inliers], 1)
        res2 = np.abs(ya_arr - (m * xa_arr + b))
        inliers2 = res2 <= 3.5
        if inliers2.sum() >= 3:
            m, b = np.polyfit(xa_arr[inliers2], ya_arr[inliers2], 1)
    return float(m * col_L + b), float(m * col_R + b), float(m)


def _find_s11_caps(mask_raw: np.ndarray, box: tuple[int, int, int, int],
                   cx: int, row_v: float, is_right_border: bool,
                   max_d: float = 540.0) -> tuple[float | None, float | None, str]:
    """Locate symmetric horizontal error-bar caps in this color's mask at column cx around row_v."""
    c0, r0, c1, r1 = box
    H, W = mask_raw.shape
    r_lo = max(r0 + 5, int(round(row_v - max_d)))
    r_hi = min(r1 + 4, int(round(row_v + max_d)))

    # MATLAB error-bar caps in S11 extend ~13 px left and right of cx; start at 5 px to clear central stems
    L_arm = mask_raw[r_lo:r_hi, max(c0 + 2, cx - 13):max(c0 + 3, cx - 5)].mean(axis=1)
    R_arm = mask_raw[r_lo:r_hi, min(W - 1, cx + 6):min(W, cx + 14)].mean(axis=1)
    cap_both_score = np.minimum(L_arm, R_arm)
    cap_any_score = np.maximum(L_arm, R_arm)

    def extract_caps(score_arr: np.ndarray, thresh: float) -> list[float]:
        rows = np.flatnonzero(score_arr >= thresh) + r_lo
        out = []
        if rows.size:
            sp = np.flatnonzero(np.diff(rows) > 2)
            st = np.concatenate(([0], sp + 1))
            en = np.concatenate((sp, [rows.size - 1]))
            for a, b in zip(st, en):
                h = int(rows[b] - rows[a] + 1)
                if 2 <= h <= 16:
                    rc = 0.5 * (rows[a] + rows[b]) + 0.5
                    if abs(rc - row_v) >= 6.0:
                        out.append(rc)
        return out

    caps_both = extract_caps(cap_both_score, 0.50)
    above = [rc for rc in caps_both if rc < row_v - 6.0]
    below = [rc for rc in caps_both if rc > row_v + 6.0]

    # Prefer symmetric pair among both-arm caps
    best_pair = None
    best_score = 1e9
    for rt in above:
        dt = row_v - rt
        for rb in below:
            db = rb - row_v
            asym = abs(dt - db)
            if asym <= max(8.0, 0.08 * 0.5 * (dt + db)):
                if asym < best_score:
                    best_score = asym
                    best_pair = (rt, rb)
    if best_pair is not None:
        return best_pair[0], best_pair[1], "both caps"

    # Match one both-arm cap with one single-arm cap (if the opposite arm is half-occluded by another color)
    caps_any = extract_caps(cap_any_score, 0.60)
    best_pair2 = None
    best_score2 = 1e9
    for rt in above:
        dt = row_v - rt
        for rb in caps_any:
            if rb > row_v + 6.0:
                asym = abs((rb - row_v) - dt)
                if asym <= max(8.0, 0.08 * dt) and asym < best_score2:
                    best_score2 = asym
                    best_pair2 = (rt, rb)
    for rb in below:
        db = rb - row_v
        for rt in caps_any:
            if rt < row_v - 6.0:
                asym = abs((row_v - rt) - db)
                if asym <= max(8.0, 0.08 * db) and asym < best_score2:
                    best_score2 = asym
                    best_pair2 = (rt, rb)
    if best_pair2 is not None:
        return best_pair2[0], best_pair2[1], "both caps (one half-occluded)"

    # Fallback: single clean both-arm cap mirrored symmetrically around row_v
    if above or below:
        stem_col = mask_raw[:, max(0, cx - 4):min(W, cx + 5)].any(axis=1)
        valid = []
        for rc in above:
            d = row_v - rc
            seg = stem_col[int(rc) + 4:min(int(row_v), int(rc) + 30)]
            if seg.size and seg.mean() >= 0.45:
                valid.append((d, rc, row_v + d, "upper cap mirrored"))
        for rc in below:
            d = rc - row_v
            seg = stem_col[max(int(row_v), int(rc) - 30):int(rc) - 4]
            if seg.size and seg.mean() >= 0.45:
                valid.append((d, row_v - d, rc, "lower cap mirrored"))
        if valid:
            valid.sort(key=lambda t: -t[0])
            _, rt, rb, tag = valid[0]
            return rt, rb, tag
        cands = [(row_v - rc, rc, row_v + (row_v - rc), "upper cap mirrored") for rc in above] + \
                [(rc - row_v, row_v - (rc - row_v), rc, "lower cap mirrored") for rc in below]
        cands.sort(key=lambda t: -t[0])
        _, rt, rb, tag = cands[0]
        return rt, rb, tag

    # Final fallback for right-border caps where left arm is occluded by another curve: use outside-right arm!
    if is_right_border:
        caps_r = extract_caps(R_arm, 0.55)
        ab_r = [rc for rc in caps_r if rc < row_v - 10.0]
        bl_r = [rc for rc in caps_r if rc > row_v + 10.0]
        for rt in ab_r:
            for rb in bl_r:
                if abs((row_v - rt) - (rb - row_v)) <= max(9.0, 0.08 * 0.5 * (rb - rt)):
                    return rt, rb, "both right-arm caps"
        if ab_r or bl_r:
            cands = [(row_v - rc, rc, row_v + (row_v - rc), "upper right-arm mirrored") for rc in ab_r] + \
                    [(rc - row_v, row_v - (rc - row_v), rc, "lower right-arm mirrored") for rc in bl_r]
            cands.sort(key=lambda t: -t[0])
            _, rt, rb, tag = cands[0]
            return rt, rb, tag

    return None, None, "no cap"


def read_s11_series(image: np.ndarray, panel: str, ma: float, at: list[float]) -> list[dict]:
    """Read one M:A curve + symmetric SD error bars from S11 panel a..f."""
    info, xa, ya = calibrate_s11(image)
    has_legend = panel in ("c", "f")
    color_map = {0.0225: "red", 0.09: "blue", 0.225: "green", 0.675: "black"}
    color = color_map[ma]

    mask_raw = _raw_color_mask(image, color, info["box"], has_legend=has_legend)
    mask_line = _line_mask(mask_raw, color, info["box"])
    c0, r0, c1, r1 = info["box"]

    x_cols = {xv: xa.to_pixel(xv) for xv in [0.1, 0.3, 0.6, 1.0, 3.0, 6.0, 10.0]}
    seg_13 = _fit_line_between(mask_line, x_cols[1.0], x_cols[3.0], pad_L=28, pad_R=28)
    seg_36 = _fit_line_between(mask_line, x_cols[3.0], x_cols[6.0], pad_L=28, pad_R=28)
    seg_610 = _fit_line_between(mask_line, x_cols[6.0], x_cols[10.0], pad_L=28, pad_R=28)

    verts: dict[float, float] = {}
    if seg_13 is not None:
        verts[1.0] = seg_13[0]
        verts[3.0] = seg_13[1]
    if seg_36 is not None:
        verts[3.0] = 0.5 * (verts[3.0] + seg_36[0]) if 3.0 in verts else seg_36[0]
        verts[6.0] = seg_36[1]
    if seg_610 is not None:
        verts[6.0] = 0.5 * (verts[6.0] + seg_610[0]) if 6.0 in verts else seg_610[0]
        verts[10.0] = seg_610[1]

    if color == "black" and 1.0 in verts:
        y_at_1 = ya.to_value(verts[1.0])
        for xv in (0.6, 0.3):
            cx = int(round(x_cols[xv]))
            r_est = ya.to_pixel(y_at_1 * xv)
            r_lo = max(r0 + 10, int(round(r_est - 100)))
            r_hi = min(r1 + 4, int(round(r_est + 100)))
            L_arm = mask_raw[r_lo:r_hi, max(c0 + 2, cx - 11):max(c0 + 3, cx - 4)].mean(axis=1)
            R_arm = mask_raw[r_lo:r_hi, min(c1 - 1, cx + 5):min(c1, cx + 12)].mean(axis=1)
            rows = np.flatnonzero(np.minimum(L_arm, R_arm) >= 0.55) + r_lo
            caps = []
            if rows.size:
                sp = np.flatnonzero(np.diff(rows) > 2)
                st = np.concatenate(([0], sp + 1))
                en = np.concatenate((sp, [rows.size - 1]))
                for a, b in zip(st, en):
                    if 2 <= (rows[b] - rows[a] + 1) <= 15:
                        caps.append(0.5 * (rows[a] + rows[b]) + 0.5)
            above = [rc for rc in caps if rc < r_est + 15]
            below = [rc for rc in caps if rc > r_est - 15]
            if above and below and min(above) < max(below) - 8:
                verts[xv] = 0.5 * (min(above) + max(below))
            elif above:
                verts[xv] = 0.5 * (min(above) + info["spines"]["bottom"])

    xs = at if at else (CHI_BLACK if color == "black" else CHI_COLOR)
    out = []
    for xv in xs:
        match_x = [k for k in verts if abs(k - xv) < 1e-4]
        if not match_x:
            out.append({"x": float(xv), "y": None, "note": f"S11({panel}) M:A={ma}: chi={xv} not found"})
            continue
        xk = match_x[0]
        rv = verts[xk]
        pc = x_cols[xk]
        cx = int(round(pc))

        max_d = 120.0 if xk <= 1.0 else (240.0 if xk <= 3.0 else (360.0 if xk <= 6.0 else 540.0))
        rt, rb, tag = _find_s11_caps(
            mask_raw, info["box"], cx, rv, is_right_border=(xk == 10.0), max_d=max_d,
        )
        if rt is not None and rb is not None and "both" in tag:
            mid = 0.5 * (rt + rb)
            if abs(mid - rv) <= 8.0:
                rv = 0.5 * (rv + mid)

        note = f"S11({panel}) M:A={ma:g} ({color}) at chi={xk:g} ({tag}; rx={xa.residual_px:.2f}px, ry={ya.residual_px:.2f}px)"
        out.append(
            C.make_reading(
                xk, rv, ya, pc, half_px=2.5, x_axis=xa, note=note,
                bar_top_px=rt, bar_bot_px=rb, bar_err_px=3.0, clamp_y=(0.0, 7.0),
            )
        )
    return out
