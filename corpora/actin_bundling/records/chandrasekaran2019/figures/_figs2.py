"""Calibrated image reader for S2 Figure (p. 23 / SI-4, mislabeled 'S1 Figure'): Shape Parameter (Sh)
vs Nematic Order Parameter (S, misspelled 'Nemetic') for the cleanly separated blue Aster-Like (AL)
cluster in (a) Unipolar and (b) Apolar non-treadmilling networks, including vertical Sh SD error bars
(900 dpi).

Dense overlapping red (ABI) and green (BL) clusters in the upper-right corner (and their uncalibrated
zoomed insets lacking numeric tick labels) are documented in figures.md.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from scipy import ndimage
import _common as C


def calibrate_s2(image: np.ndarray):
    info = C.find_matlab_box(
        image, dark_thresh=140.0, min_h_frac=0.55, min_v_frac=0.55,
        has_top=True, has_right=True, min_tick_len=5, max_tick_len=50,
        max_tick_w=30,
    )
    sp = info["spines"]
    xt = C.merge_tick_set(sp["left"], info["x_ticks_top"] + info["x_ticks_bot"], sp["right"], 6, tol_px=25.0)
    yt = C.merge_tick_set(sp["top"], info["y_ticks_left"] + info["y_ticks_right"], sp["bottom"], 6, tol_px=25.0)
    xa = C.fit_axis(xt, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], max_residual_px=4.0)
    ya = C.fit_axis(yt, [1.0, 0.8, 0.6, 0.4, 0.2, 0.0], max_residual_px=4.0)
    return info, xa, ya


def _find_vertical_sd_bar(dark: np.ndarray, blue_mask: np.ndarray, cx: int, row_c: float,
                          r_lo: int, r_hi: int, marker_rad: int = 24,
                          max_d: int = 260) -> tuple[float | None, float | None, str]:
    """Trace the black vertical SD error bar extending above and below a blue AL circle at (cx, row_c)."""
    H, W = dark.shape
    stem = dark[:, max(0, cx - 4):min(W, cx + 5)].any(axis=1)
    other_blue = blue_mask[:, max(0, cx - 18):min(W, cx + 19)].any(axis=1)

    y0 = max(r_lo, int(round(row_c - max_d)))
    y1 = min(r_hi, int(round(row_c + max_d)))
    L_arm = dark[y0:y1, max(0, cx - 12):max(1, cx - 5)].mean(axis=1)
    R_arm = dark[y0:y1, min(W - 1, cx + 6):min(W, cx + 13)].mean(axis=1)
    cap_rows = np.flatnonzero(np.minimum(L_arm, R_arm) >= 0.45) + y0
    caps = []
    if cap_rows.size:
        sp = np.flatnonzero(np.diff(cap_rows) > 2)
        st = np.concatenate(([0], sp + 1))
        en = np.concatenate((sp, [cap_rows.size - 1]))
        for a, b in zip(st, en):
            h = int(cap_rows[b] - cap_rows[a] + 1)
            if 2 <= h <= 14:
                rc = 0.5 * (cap_rows[a] + cap_rows[b]) + 0.5
                if abs(rc - row_c) >= marker_rad + 10:
                    caps.append(rc)

    def has_clean_stem(rc: float) -> bool:
        a, b = int(round(min(row_c, rc))), int(round(max(row_c, rc)))
        if b - a <= marker_rad:
            return False
        if rc < row_c:
            seg_s = stem[a:int(round(row_c)) - marker_rad]
            seg_b = other_blue[a:int(round(row_c)) - marker_rad]
        else:
            seg_s = stem[int(round(row_c)) + marker_rad:b + 1]
            seg_b = other_blue[int(round(row_c)) + marker_rad:b + 1]
        if seg_b.sum() > 2:
            return False
        return bool(seg_s.size > 0 and seg_s.mean() >= 0.85)

    above = [rc for rc in caps if rc < row_c - marker_rad - 10 and has_clean_stem(rc)]
    below = [rc for rc in caps if rc > row_c + marker_rad + 10 and has_clean_stem(rc)]

    best = None
    best_asym = 1e9
    for rt in above:
        dt = row_c - rt
        for rb in below:
            db = rb - row_c
            asym = abs(dt - db)
            if asym <= max(10.0, 0.12 * 0.5 * (dt + db)) and asym < best_asym:
                best_asym = asym
                best = (rt, rb)
    if best is not None:
        return best[0], best[1], "both caps"

    single = []
    for rt in above:
        d = row_c - rt
        single.append((d, rt, row_c + d, "upper cap mirrored"))
    for rb in below:
        d = rb - row_c
        single.append((d, row_c - d, rb, "lower cap mirrored"))
    if single:
        # Prefer the cap whose stem terminates past it, else closest clean cap
        single.sort(key=lambda t: (
            0 if (
                (t[3].startswith("upper") and stem[max(r_lo, int(t[1]) - 18):max(r_lo, int(t[1]) - 6)].mean() < 0.25)
                or (t[3].startswith("lower") and stem[min(r_hi, int(t[2]) + 7):min(r_hi, int(t[2]) + 19)].mean() < 0.25)
            ) else 1,
            t[0],
        ))
        _, rt, rb, tag = single[0]
        return rt, rb, tag
    return None, None, "no cap"


def read_s2_al(image: np.ndarray, polarity: str, at: list[float]) -> list[dict]:
    """Read cleanly resolved blue Aster-Like (AL) points + vertical Sh SD bars from S2 Figure."""
    info, xa, ya = calibrate_s2(image)
    c0, r0, c1, r1 = info["box"]

    r = image[:, :, 0].astype(np.int16)
    g = image[:, :, 1].astype(np.int16)
    b = image[:, :, 2].astype(np.int16)
    dark = (r + g + b) < 330

    mask_main = np.zeros(image.shape[:2], dtype=bool)
    mask_main[r0 + 8:r1 - 8, c0 + 8:c1 - 8] = True
    mask_main[int(ya.to_pixel(0.51)):r1, int(xa.to_pixel(0.55)):c1] = False

    m_blue = mask_main & (b > 160) & (r < 90) & (g < 100) & ((b - r) > 90)
    lab, n = ndimage.label(m_blue)

    min_area = 1200
    pts = []
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        if rr.size < min_area:
            continue
        rc = float(rr.mean()) + 0.5
        cc_m = float(cc.mean()) + 0.5
        s_val = xa.to_value(cc_m)
        cx = int(round(cc_m))
        r_hi_eff = int(ya.to_pixel(0.505)) if s_val > 0.55 else (r1 - 5)
        rt, rb, tag = _find_vertical_sd_bar(dark, m_blue & (lab != i), cx, rc, r0 + 5, r_hi_eff)
        pts.append((s_val, rc, cc_m, int(rr.size), rt, rb, tag))

    pts.sort(key=lambda t: t[0])
    out = []
    for s_val, rc, cc_m, area, rt, rb, tag in pts:
        note = (
            f"S2 {polarity} AL blue circle (area={area}px, {tag}; "
            f"rx={xa.residual_px:.2f}px, ry={ya.residual_px:.2f}px)"
        )
        out.append(
            C.make_reading(
                s_val, rc, ya, cc_m, half_px=4.0, x_axis=xa, note=note,
                bar_top_px=rt, bar_bot_px=rb, bar_err_px=4.0, clamp_y=(0.0, 1.0),
            )
        )
    return out
