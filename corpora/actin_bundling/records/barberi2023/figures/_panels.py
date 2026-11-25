"""Panel-family figure readers for barberi2023.

Each function calibrates directly from the rendered crop `img` (RGB numpy array),
extracts the target curve or markers, computes honest 95% intervals including
line thickness, calibration misfit, and slope-propagated horizontal error, and
returns BoxGraph `mode: code` reading dicts with pixel coordinates `px`.
"""

from __future__ import annotations

import math
import numpy as np
from scipy import ndimage

import _common as C

L_HALF = 5.0 * math.pi  # L/2 = 10*pi / 2 = 15.707963267948966


def _filter_components_by_span(mask: np.ndarray, lo: int, hi: int | None = None) -> np.ndarray:
    """Keep connected components whose max(height, width) lies in [lo, hi]."""
    lab, n = ndimage.label(mask)
    if n == 0:
        return np.zeros_like(mask, dtype=bool)
    objs = ndimage.find_objects(lab)
    out = np.zeros_like(mask, dtype=bool)
    for idx, sl in enumerate(objs, start=1):
        if sl is None:
            continue
        span = max(sl[0].stop - sl[0].start, sl[1].stop - sl[1].start)
        if span >= lo and (hi is None or span <= hi):
            out[sl] |= (lab[sl] == idx)
    return out


def _mask_rel_rect(mask: np.ndarray, frame: dict, x0f: float, x1f: float, y0f: float, y1f: float) -> None:
    """Zero out a relative rectangle inside the plot box (for panel labels / legends)."""
    c0, r0, c1, r1 = frame["box"]
    W = c1 - c0
    H = r1 - r0
    ca = max(0, int(round(c0 + x0f * W)))
    cb = min(mask.shape[1], int(round(c0 + x1f * W)))
    ra = max(0, int(round(r0 + y0f * H)))
    rb = min(mask.shape[0], int(round(r0 + y1f * H)))
    mask[ra:rb, ca:cb] = False


# ==============================================================================
# 1. Main Text Fig. 1(a): LP0 (C, Na, Ni, V)
# ==============================================================================
def read_fig1a(img: np.ndarray, at: list[float], curve: str) -> list[dict]:
    fr = C.find_frame(img, include_pink_right=True)
    bt = C.find_bottom_ticks(img, fr)
    lt = C.find_left_ticks(img, fr)
    rt = C.find_right_ticks(img, fr, use_pink=True)
    if len(bt) < 2 or len(lt) < 2 or len(rt) < 2:
        raise ValueError(f"Fig. 1(a) ticks missing: bt={bt}, lt={lt}, rt={rt}")
    xa = C.LinearAxis([(bt[0], -L_HALF), (bt[-1], L_HALF)])
    ya_left = C.LinearAxis([(lt[0], 3.0), (lt[-1], 0.0)])
    ya_right = C.LinearAxis([(rt[0], 0.7), (rt[-1], -0.7)])

    if curve == "V":
        m = C.interior_mask(C.pink_mask(img), fr, inset=6)
        return C.sample_curve_column(m, xa, ya_right, at, half_win=2, gap_search=8, note_prefix="Fig.1a V: ")

    dm = C.interior_mask(C.dark_mask(img, 145.0), fr, inset=6)
    # Remove '(a)' top-left, legend top-right, 'LP0' label (strictly above C baseline at y~0.945)
    _mask_rel_rect(dm, fr, 0.0, 0.24, 0.0, 0.11)
    _mask_rel_rect(dm, fr, 0.56, 1.0, 0.0, 0.20)
    _mask_rel_rect(dm, fr, 0.76, 0.97, 0.82, 0.915)

    m_solid = _filter_components_by_span(dm, lo=160)
    if curve == "C":
        return C.sample_curve_column(m_solid, xa, ya_left, at, half_win=2, gap_search=14, note_prefix="Fig.1a C: ")
    if curve == "Na":
        m_dash = _filter_components_by_span(dm, lo=15, hi=55)
        res = C.sample_curve_column(m_dash, xa, ya_left, at, half_win=3, gap_search=40, note_prefix="Fig.1a Na: ")
        # In outer flat wings (|X| >= 4), Na == C identically so dashed Na merges with solid C baseline
        for i, r in enumerate(res):
            if r.get("y") is None:
                res[i] = C.sample_curve_column(
                    m_solid, xa, ya_left, [at[i]], half_win=2, gap_search=14, note_prefix="Fig.1a Na(wing=C): "
                )[0]
        return res
    if curve == "Ni":
        m = _filter_components_by_span(dm, lo=3, hi=14)
        return C.sample_curve_column(m, xa, ya_left, at, half_win=3, gap_search=28, branch="top", note_prefix="Fig.1a Ni: ")
    raise ValueError(f"Unknown Fig. 1(a) curve: {curve}")


# ==============================================================================
# 2. Main Text Fig. 1(b)-(g): Multipeaked LP2..LP7 actomyosin C(X)
# ==============================================================================
def read_fig1_lp(img: np.ndarray, at: list[float], label: str = "LP") -> list[dict]:
    fr = C.find_frame(img, include_pink_right=False)
    bt = C.find_bottom_ticks(img, fr, pad_lo=6, pad_hi=26)
    lt = C.find_left_ticks(img, fr, pad_lo=6, pad_hi=26)
    if len(bt) < 2 or len(lt) < 2:
        raise ValueError(f"{label} ticks missing: bt={bt}, lt={lt}")
    xa = C.LinearAxis([(bt[0], -L_HALF), (bt[-1], L_HALF)])
    ya = C.LinearAxis([(lt[0], 1.0), (lt[-1], 0.0)])

    dm = C.interior_mask(C.dark_mask(img, 145.0), fr, inset=6)
    m = C.largest_span_component(dm, min_span_frac=0.55)
    return C.sample_curve_column(m, xa, ya, at, half_win=2, gap_search=4, note_prefix=f"{label} C: ")


# ==============================================================================
# 3. Main Text Fig. 2: Slanted Snaking Bifurcation Diagram
# ==============================================================================
def _detect_dots_edt(mask: np.ndarray, min_r: float = 12.0, foot: int = 15, merge_px: float = 16.0) -> list[tuple[float, float, float]]:
    """Fast circular marker center detector using 2x downsampling + separable 1D max filter."""
    # Downsample by 2 for 16x speedup, then refine centroids on full-res local windows
    sub = mask[::2, ::2]
    dist = ndimage.distance_transform_edt(sub)
    f2 = max(3, foot // 2 | 1)
    mx = ndimage.maximum_filter1d(ndimage.maximum_filter1d(dist, size=f2, axis=0), size=f2, axis=1)
    loc = (dist == mx) & (dist >= min_r * 0.5)
    lab, n = ndimage.label(loc)
    objs = ndimage.find_objects(lab)
    raw = []
    H, W = mask.shape
    for idx, sl in enumerate(objs, start=1):
        if sl is None:
            continue
        rr, cc = np.nonzero(lab[sl] == idx)
        r0 = (float(rr.mean()) + sl[0].start) * 2.0 + 1.0
        c0 = (float(cc.mean()) + sl[1].start) * 2.0 + 1.0
        rad = float(dist[int(rr.mean()) + sl[0].start, int(cc.mean()) + sl[1].start]) * 2.0
        # Sub-pixel centroid refinement inside inscribed disk of radius 0.75 * rad
        r_ref = max(4, int(round(0.7 * rad)))
        ra, rb = max(0, int(r0) - r_ref), min(H, int(r0) + r_ref + 1)
        ca, cb = max(0, int(c0) - r_ref), min(W, int(c0) + r_ref + 1)
        patch = mask[ra:rb, ca:cb]
        pr, pc = np.nonzero(patch)
        if len(pr) > 10:
            r0 = float(pr.mean()) + ra + 0.5
            c0 = float(pc.mean()) + ca + 0.5
        raw.append((c0, r0, rad))
    raw.sort(key=lambda t: -t[2])
    kept: list[tuple[float, float, float]] = []
    for c, r, rad in raw:
        if all(math.hypot(c - kc, r - kr) >= merge_px for kc, kr, _ in kept):
            kept.append((c, r, rad))
    return kept


def read_fig2_branch(img: np.ndarray, at: list[float], branch: str) -> list[dict]:
    fr = C.find_frame(img, include_pink_right=False)
    bt = C.find_bottom_ticks(img, fr)
    lt = C.find_left_ticks(img, fr)
    if len(bt) < 2 or len(lt) < 2:
        raise ValueError(f"Fig. 2 ticks missing: bt={bt}, lt={lt}")
    xa = C.LinearAxis([(bt[0], 6.0), (bt[-1], 19.0)])
    ya = C.LinearAxis([(lt[0], 0.12), (lt[-1], 0.0)])

    R = img[:, :, 0].astype(np.int16)
    G = img[:, :, 1].astype(np.int16)
    B = img[:, :, 2].astype(np.int16)

    if branch == "H":
        m = C.interior_mask((R > 140) & (G < 80) & (B > 40) & (R - G > 80), fr, inset=4)
        pts = _detect_dots_edt(m, min_r=14.0, foot=17, merge_px=20.0)
        selected = sorted([(xa.to_val(c), ya.to_val(r), c, r, rad) for c, r, rad in pts])
    elif branch == "LP0":
        m = C.interior_mask((G > R + 35) & (G > B + 35) & (G > 60), fr, inset=4)
        c_cut = int(round(xa.to_px(8.0)))
        m[:, c_cut:] = False
        if at:
            return C.sample_curve_column(m, xa, ya, at, half_win=2, gap_search=4, min_half_px=2.0, note_prefix="Fig.2 LP0: ")
        pts = _detect_dots_edt(m, min_r=14.0, foot=9, merge_px=12.0)
        selected = sorted([(xa.to_val(c), ya.to_val(r), c, r, rad) for c, r, rad in pts])
    elif branch in ("LP2", "LP4", "LP6"):
        m = C.interior_mask((B > R + 60) & (B > G + 40) & (B > 100), fr, inset=4)
        pts = _detect_dots_edt(m, min_r=14.0, foot=17, merge_px=20.0)
        vals = sorted([(xa.to_val(c), ya.to_val(r), c, r, rad) for c, r, rad in pts])
        if branch == "LP2":
            selected = [p for p in vals if p[0] <= 10.5]
        elif branch == "LP4":
            selected = [p for p in vals if 10.6 <= p[0] <= 14.2 and (p[0] < 13.2 or p[1] > 0.075)]
        else:  # LP6
            selected = [p for p in vals if p[0] >= 13.2 and (p[0] > 14.2 or p[1] < 0.075)]
    elif branch in ("LP3", "LP5", "LP7"):
        m = C.interior_mask((G > R + 35) & (G > B + 35) & (G > 60), fr, inset=4)
        pts = _detect_dots_edt(m, min_r=14.0, foot=17, merge_px=20.0)
        vals = sorted([(xa.to_val(c), ya.to_val(r), c, r, rad) for c, r, rad in pts])
        if branch == "LP3":
            selected = [p for p in vals if 9.0 <= p[0] <= 12.5]
        elif branch == "LP5":
            selected = [p for p in vals if 12.6 <= p[0] <= 17.2 and (p[0] < 15.2 or p[1] > 0.085)]
        else:  # LP7
            selected = [p for p in vals if p[0] >= 15.2 and (p[0] > 17.2 or p[1] < 0.085)]
    else:
        raise ValueError(f"Unknown Fig. 2 branch: {branch}")

    y_base_px = 3.0 + ya.residual_px
    x_err_px = 2.0 + xa.residual_px
    out = []
    for i, (x, y, c, r, rad) in enumerate(selected):
        # Local branch slope |dr/dc| from neighboring dots so horizontal calibration uncertainty propagates honestly
        slope = 0.0
        if len(selected) >= 2:
            if 0 < i < len(selected) - 1:
                dc = abs(selected[i + 1][2] - selected[i - 1][2])
                dr = abs(selected[i + 1][3] - selected[i - 1][3])
                slope = dr / max(1.0, dc)
            elif i == 0:
                slope = abs(selected[1][3] - r) / max(1.0, abs(selected[1][2] - c))
            else:
                slope = abs(r - selected[-2][3]) / max(1.0, abs(c - selected[-2][2]))
        half_px = y_base_px + slope * x_err_px
        ya_a, ya_b = ya.to_val(r - half_px), ya.to_val(r + half_px)
        lo, hi = min(ya_a, ya_b), max(ya_a, ya_b)
        out.append({
            "x": float(x),
            "y": float(y),
            "lo": float(lo),
            "hi": float(hi),
            "px": [float(c), float(r)],
            "px_lo": float(r - half_px),
            "px_hi": float(r + half_px),
            "note": f"Fig.2 {branch} dot (r={rad:.1f}px, slope={slope:.2f})",
        })
    return out


# ==============================================================================
# 4. Main Text Fig. 3(a): Spike profile (C, Na, Ni, V)
# ==============================================================================
def read_fig3a(img: np.ndarray, at: list[float], curve: str) -> list[dict]:
    fr = C.find_frame(img, include_pink_right=True)
    bt = C.find_bottom_ticks(img, fr)
    lt = C.find_left_ticks(img, fr)
    rt = C.find_right_ticks(img, fr, use_pink=True)
    if len(bt) < 2 or len(lt) < 2 or len(rt) < 2:
        raise ValueError(f"Fig. 3(a) ticks missing: bt={bt}, lt={lt}, rt={rt}")
    xa = C.LinearAxis([(bt[0], -L_HALF), (bt[-1], L_HALF)])
    ya_left = C.LinearAxis([(lt[0], 1.2), (lt[-1], 0.0)])
    ya_right = C.LinearAxis([(rt[0], 0.2), (rt[-1], -0.2)])

    if curve == "V":
        m = C.interior_mask(C.pink_mask(img), fr, inset=6)
        return C.sample_curve_column(m, xa, ya_right, at, half_win=2, gap_search=8, note_prefix="Fig.3a V: ")

    dm = C.interior_mask(C.dark_mask(img, 145.0), fr, inset=6)
    _mask_rel_rect(dm, fr, 0.0, 0.22, 0.0, 0.20)   # '(a)'
    _mask_rel_rect(dm, fr, 0.56, 0.98, 0.48, 0.85)  # legend C, Na, Ni

    if curve == "C":
        m = _filter_components_by_span(dm, lo=160)
        return C.sample_curve_column(m, xa, ya_left, at, half_win=1, gap_search=10, branch="top", note_prefix="Fig.3a C: ")
    if curve == "Na":
        m_dash = _filter_components_by_span(dm, lo=15, hi=55)
        m_solid = _filter_components_by_span(dm, lo=160)
        res = C.sample_curve_column(m_dash, xa, ya_left, at, half_win=3, gap_search=32, note_prefix="Fig.3a Na: ")
        for i, r in enumerate(res):
            if r.get("y") is None and abs(at[i]) > 1.8:
                res[i] = C.sample_curve_column(
                    m_solid, xa, ya_left, [at[i]], half_win=2, gap_search=14, branch="bot", note_prefix="Fig.3a Na(wing=C): "
                )[0]
        return res
    if curve == "Ni":
        m = _filter_components_by_span(dm, lo=3, hi=14)
        return C.sample_curve_column(m, xa, ya_left, at, half_win=3, gap_search=28, branch="top", note_prefix="Fig.3a Ni: ")
    raise ValueError(f"Unknown Fig. 3(a) curve: {curve}")


# ==============================================================================
# 5. Main Text Fig. 3(b): Spatial Dynamics Phase Portrait C vs d_X C
# ==============================================================================
def read_fig3b(img: np.ndarray, at: list[float], orbit: str, half: str) -> list[dict]:
    fr = C.find_frame(img, include_pink_right=False)
    bt = C.find_bottom_ticks(img, fr)
    lt = C.find_left_ticks(img, fr)
    if len(bt) < 2 or len(lt) < 2:
        raise ValueError(f"Fig. 3(b) ticks missing: bt={bt}, lt={lt}")
    xa = C.LinearAxis([(bt[0], -4.0), (bt[-1], 4.0)])
    ya = C.LinearAxis([(lt[0], 0.8), (lt[-1], 0.2)])

    dm = C.interior_mask(C.dark_mask(img, 145.0), fr, inset=6)
    _mask_rel_rect(dm, fr, 0.0, 0.22, 0.0, 0.20)  # '(b)'

    if orbit == "spike":
        m = _filter_components_by_span(dm, lo=250)
        br = "top" if half == "upper" else "bot"
        return C.sample_curve_column(m, xa, ya, at, half_win=2, gap_search=6, branch=br, note_prefix=f"Fig.3b spike {half}: ")
    if orbit == "LP0":
        # Separate upper arc (C >= 0.58) and lower arc (C <= 0.58) so gaps near the crossing never jump halves
        m_dash = _filter_components_by_span(dm, lo=15, hi=60)
        r_split = int(round(ya.to_px(0.58)))
        if half == "upper":
            # Above C = 0.68, LP0 is completely above the solid spike curve, so full dm also recovers crossing dashes
            m = m_dash.copy()
            r_above_spike = int(round(ya.to_px(0.665)))
            m[:r_above_spike, :] = dm[:r_above_spike, :]
            m[r_split:, :] = False
            return C.sample_curve_column(m, xa, ya, at, half_win=3, gap_search=40, branch="top", note_prefix="Fig.3b LP0 upper: ")
        else:
            m = m_dash.copy()
            m[:r_split, :] = False
            return C.sample_curve_column(m, xa, ya, at, half_win=3, gap_search=40, branch="bot", note_prefix="Fig.3b LP0 lower: ")
    raise ValueError(f"Unknown Fig. 3(b) orbit: {orbit}")


# ==============================================================================
# 6. Main Text Fig. 3(c): Broken-Axis Spatial Eigenvalues
# ==============================================================================
def read_fig3c(img: np.ndarray, at: list[float], which: str) -> list[dict]:
    """Read plotted x (spike) and + (LP0) eigenvalues from Fig. 3(c) broken axes."""
    fr = C.find_frame(img, include_pink_right=False, min_frac=0.35)
    lt = C.find_left_ticks(img, fr)
    bt = C.find_bottom_ticks(img, fr)
    if len(lt) != 5:
        raise ValueError(f"Fig. 3(c) expected 5 y ticks (20,10,0,-10,-20), got {lt}")
    ya = C.LinearAxis(list(zip(lt, [20.0, 10.0, 0.0, -10.0, -20.0])))

    if len(bt) < 6:
        raise ValueError(f"Fig. 3(c) bottom ticks unexpected: {bt}")
    xa_left = C.LinearAxis([(bt[0], 0.0), (bt[1], 1.0)])
    xa_right = C.LinearAxis(list(zip(bt[-4:], [40.0, 60.0, 80.0, 100.0])))
    if xa_right.residual_px > 1.5:
        raise ValueError(f"Fig. 3(c) right x-axis misfit {xa_right.residual_px:.2f} px on {bt[-4:]}")

    dark = C.dark_mask(img, thresh=95.0)
    c0, r0, c1, r1 = fr["box"]
    c_mid_l = int(round(bt[2] - 20)) if len(bt) >= 8 else int(round(0.45 * (c0 + c1)))
    c_mid_r = int(round(bt[-5] + 20)) if len(bt) >= 8 else int(round(0.55 * (c0 + c1)))

    left_box = np.zeros_like(dark)
    left_box[r0 + 8 : r1 - 8, c0 + 10 : c_mid_l] = dark[r0 + 8 : r1 - 8, c0 + 10 : c_mid_l]
    left_box[r0 : r0 + int(0.18 * (r1 - r0)), c0 : c0 + int(0.10 * (c1 - c0))] = False

    right_box = np.zeros_like(dark)
    right_box[r0 + 8 : r1 - 8, c_mid_r : c1 - 10] = dark[r0 + 8 : r1 - 8, c_mid_r : c1 - 10]

    def _find_markers(bmask, xa_side, side_name):
        lab, n = ndimage.label(ndimage.binary_dilation(bmask, iterations=2))
        objs = ndimage.find_objects(lab)
        pts = []
        for idx, sl in enumerate(objs, start=1):
            if sl is None:
                continue
            rr, cc = np.nonzero((lab[sl] == idx) & bmask[sl])
            if len(rr) < 400:
                continue
            r = float(rr.mean()) + sl[0].start + 0.5
            c = float(cc.mean()) + sl[1].start + 0.5
            re_v = xa_side.to_val(c)
            im_v = ya.to_val(r)
            ir, ic = int(round(r)), int(round(c))
            h_arm = int(bmask[max(0, ir - 2) : ir + 3, max(0, ic - 25) : ic + 26].sum())
            v_arm = int(bmask[max(0, ir - 25) : ir + 26, max(0, ic - 2) : ic + 3].sum())
            is_plus = (h_arm > 180 and v_arm > 180)
            is_both = (len(rr) > 2600)
            pts.append({
                "re": re_v, "im": im_v, "c": c, "r": r,
                "re_err": xa_side.val_per_px * (2.5 + xa_side.residual_px),
                "im_err": ya.val_per_px * (2.5 + ya.residual_px),
                "is_plus": is_plus, "is_both": is_both, "side": side_name,
            })
        return pts

    all_pts = _find_markers(left_box, xa_left, "L") + _find_markers(right_box, xa_right, "R")
    spike_pts = sorted([p for p in all_pts if p["is_both"] or not p["is_plus"]], key=lambda p: (p["re"], p["im"]))
    lp0_pts = sorted([p for p in all_pts if p["is_both"] or p["is_plus"]], key=lambda p: (p["re"], p["im"]))

    if which == "spike_re":
        return [{
            "x": float(i + 1), "y": p["re"], "lo": p["re"] - p["re_err"], "hi": p["re"] + p["re_err"],
            "px": [p["c"], p["r"]],
            "note": f"Fig.3c spike marker {i+1}: Re={p['re']:.3f}, Im={p['im']:.3f}",
        } for i, p in enumerate(spike_pts)]
    if which == "LP0_re":
        return [{
            "x": float(i + 1), "y": p["re"], "lo": p["re"] - p["re_err"], "hi": p["re"] + p["re_err"],
            "px": [p["c"], p["r"]],
            "note": f"Fig.3c LP0 marker {i+1}: Re={p['re']:.3f}, Im={p['im']:.3f}",
        } for i, p in enumerate(lp0_pts)]
    if which == "LP0_im":
        return [{
            "x": float(i + 1), "y": p["im"], "lo": p["im"] - p["im_err"], "hi": p["im"] + p["im_err"],
            "px": [p["c"], p["r"]],
            "px_lo": p["r"] - 2.5, "px_hi": p["r"] + 2.5,
            "note": f"Fig.3c LP0 marker {i+1}: Re={p['re']:.3f}, Im={p['im']:.3f}",
        } for i, p in enumerate(lp0_pts)]
    if which == "spike_im":
        return [{
            "x": float(i + 1), "y": p["im"], "lo": p["im"] - p["im_err"], "hi": p["im"] + p["im_err"],
            "px": [p["c"], p["r"]],
            "px_lo": p["r"] - 2.5, "px_hi": p["r"] + 2.5,
            "note": f"Fig.3c spike marker {i+1}: Re={p['re']:.3f}, Im={p['im']:.3f}",
        } for i, p in enumerate(spike_pts)]
    raise ValueError(f"Unknown Fig. 3(c) mode: {which}")


# ==============================================================================
# 7. SM Fig. 1: Stability Diagram Omega_c(Z) and LS dots
# ==============================================================================
def read_smfig1(img: np.ndarray, at: list[float], kind: str) -> list[dict]:
    fr = C.find_frame(img, include_pink_right=False)
    bt = C.find_bottom_ticks(img, fr, pad_lo=6, pad_hi=22)
    lt = C.find_left_ticks(img, fr, pad_lo=6, pad_hi=22)
    if len(bt) < 2 or len(lt) < 2:
        raise ValueError(f"SM Fig. 1 ticks missing: bt={bt}, lt={lt}")
    xa = C.LinearAxis([(bt[0], 0.0), (bt[-1], 30.0)])
    ya = C.LinearAxis([(lt[0], 30.0), (lt[-1], 0.0)])

    dm = C.interior_mask(C.dark_mask(img, 120.0), fr, inset=8)
    _mask_rel_rect(dm, fr, 0.0, 0.36, 0.85, 1.0)

    # Fast dot detection via 2x downsampled EDT + separable 1D max filter (~0.15s)
    dots = _detect_dots_edt(dm, min_r=15.0, foot=31, merge_px=30.0)
    if kind == "LS_onset":
        by_z: dict[int, tuple[float, float, float, float]] = {}
        for c, r, rad in dots:
            zv = xa.to_val(c)
            omv = ya.to_val(r)
            zk = int(round(zv))
            if zk not in by_z or omv < by_z[zk][1]:
                by_z[zk] = (zv, omv, c, r)
        half_px = 3.0 + ya.residual_px
        out = []
        for zk in sorted(by_z):
            zv, omv, c, r = by_z[zk]
            a, b = ya.to_val(r - half_px), ya.to_val(r + half_px)
            out.append({
                "x": float(zv), "y": float(omv), "lo": min(a, b), "hi": max(a, b),
                "px": [c, r], "px_lo": r - half_px, "px_hi": r + half_px,
                "note": f"SM Fig.1 LS onset dot at Z~{zk}",
            })
        return out

    curve_mask = dm.copy()
    for c, r, rad in dots:
        rad_erase = rad + 10.0
        r0 = max(0, int(r - rad_erase - 1))
        r1 = min(dm.shape[0], int(r + rad_erase + 2))
        c0 = max(0, int(c - rad_erase - 1))
        c1 = min(dm.shape[1], int(c + rad_erase + 2))
        sub_y = np.arange(r0, r1)[:, None] - r
        sub_x = np.arange(c0, c1)[None, :] - c
        curve_mask[r0:r1, c0:c1][sub_x * sub_x + sub_y * sub_y <= rad_erase * rad_erase] = False

    if kind == "boundary_Io":
        return C.sample_curve_column(curve_mask, xa, ya, at, half_win=3, gap_search=30, min_half_px=2.5,
                                     note_prefix="SM Fig.1 Io: ")
    if kind == "boundary_Is":
        return C.sample_curve_column(curve_mask, xa, ya, at, half_win=3, gap_search=20, min_half_px=2.5,
                                     note_prefix="SM Fig.1 Is: ")
    raise ValueError(f"Unknown SM Fig. 1 kind: {kind}")


# ==============================================================================
# 8. SM Fig. 3 Main: Spike profile C(X) and H line over X in [-L/8, L/8]
# ==============================================================================
def read_smfig3_main(img: np.ndarray, at: list[float], kind: str) -> list[dict]:
    fr = C.find_frame(img, include_pink_right=False)
    bt = C.find_bottom_ticks(img, fr, pad_lo=5, pad_hi=22)
    lt = C.find_left_ticks(img, fr, pad_lo=5, pad_hi=22)
    if len(bt) < 2 or len(lt) < 2:
        raise ValueError(f"SM Fig. 3 main ticks missing: bt={bt}, lt={lt}")
    x_edge = 10.0 * math.pi / 8.0
    xa = C.LinearAxis([(bt[0], -x_edge), (bt[-1], x_edge)])
    ya = C.LinearAxis([(lt[0], 0.7), (lt[-1], 0.2)])

    if kind == "H":
        m = C.interior_mask(C.pink_mask(img), fr, inset=6)
        _mask_rel_rect(m, fr, 0.0, 0.20, 0.0, 0.28)
        return C.sample_curve_column(m, xa, ya, at, half_win=3, gap_search=35, note_prefix="SM Fig.3 H: ")

    if kind == "C_L32":
        gray = img.astype(np.float32).mean(axis=2)
        pink = C.pink_mask(img)
        m = C.interior_mask((gray < 175.0) & ~pink, fr, inset=6)
        _mask_rel_rect(m, fr, 0.0, 0.22, 0.0, 0.30)
        _mask_rel_rect(m, fr, 0.55, 1.0, 0.0, 0.48)
        return C.sample_curve_column(m, xa, ya, at, half_win=2, gap_search=6, branch="top", note_prefix="SM Fig.3 C(L32): ")
    raise ValueError(f"Unknown SM Fig. 3 main kind: {kind}")


# ==============================================================================
# 9. SM Fig. 3 Inset: Max. C vs System size (x L) = 1, 2, 4, 8, 16, 32
# ==============================================================================
def read_smfig3_inset(img: np.ndarray, at: list[float]) -> list[dict]:
    fr = C.find_frame(img, include_pink_right=False)
    bt = C.find_bottom_ticks(img, fr, pad_lo=6, pad_hi=28)
    lt = C.find_left_ticks(img, fr, pad_lo=6, pad_hi=28)
    if len(lt) != 4 or len(bt) != 16:
        raise ValueError(f"SM Fig. 3 inset expected 4 y-ticks & 16 x-ticks, got y={len(lt)}, x={len(bt)}")
    xa = C.LinearAxis([(p, 2.0 * (i + 1)) for i, p in enumerate(bt)])
    ya = C.LinearAxis(list(zip(lt, [0.68, 0.67, 0.66, 0.65])))
    if xa.residual_px > 1.0 or ya.residual_px > 1.0:
        raise ValueError(f"SM Fig. 3 inset misfit too large: x={xa.residual_px:.2f}, y={ya.residual_px:.2f}")

    dm = C.interior_mask(C.dark_mask(img, 100.0), fr, inset=10)
    lab, n = ndimage.label(dm)
    objs = ndimage.find_objects(lab)
    pts = []
    half_px = 2.0 + ya.residual_px
    for idx, sl in enumerate(objs, start=1):
        if sl is None:
            continue
        rr, cc = np.nonzero(lab[sl] == idx)
        if len(rr) < 1000:
            continue
        r = float(rr.mean()) + sl[0].start + 0.5
        c = float(cc.mean()) + sl[1].start + 0.5
        xv = xa.to_val(c)
        yv = ya.to_val(r)
        a, b = ya.to_val(r - half_px), ya.to_val(r + half_px)
        pts.append({
            "x": float(round(xv)) if abs(xv - round(xv)) < 0.03 else float(xv),
            "y": float(yv),
            "lo": min(a, b),
            "hi": max(a, b),
            "px": [c, r],
            "px_lo": r - half_px,
            "px_hi": r + half_px,
            "note": f"SM Fig.3 inset dot (raw x={xv:.3f}, area={len(rr)}px, misfit={ya.residual_px:.2f}px)",
        })
    pts.sort(key=lambda d: d["x"])
    return pts


# ==============================================================================
# 10. SM Fig. 4 Left (LP0) & Right (LP3): -d_X J_c, R_c, V
# ==============================================================================
def read_smfig4(img: np.ndarray, at: list[float], curve: str, label: str = "LP0") -> list[dict]:
    fr = C.find_frame(img, include_pink_right=True)
    bt = C.find_bottom_ticks(img, fr, pad_lo=5, pad_hi=22)
    lt = C.find_left_ticks(img, fr, pad_lo=5, pad_hi=22)
    rt = C.find_right_ticks(img, fr, pad_lo=5, pad_hi=22, use_pink=True)
    if len(bt) < 2 or len(lt) < 2 or len(rt) < 2:
        raise ValueError(f"SM Fig. 4 ({label}) ticks missing: bt={bt}, lt={lt}, rt={rt}")
    xa = C.LinearAxis([(bt[0], -L_HALF), (bt[-1], L_HALF)])
    ya_left = C.LinearAxis([(lt[0], 1.0), (lt[-1], -1.0)])
    ya_right = C.LinearAxis([(rt[0], 0.7), (rt[-1], -0.7)])

    if curve == "V":
        m = C.interior_mask(C.pink_mask(img), fr, inset=6)
        return C.sample_curve_column(m, xa, ya_right, at, half_win=2, gap_search=8, note_prefix=f"SM Fig.4 {label} V: ")

    dm = C.interior_mask(C.dark_mask(img, 145.0), fr, inset=6)
    _mask_rel_rect(dm, fr, 0.52, 1.0, 0.0, 0.14)  # legend top-right

    m_solid = _filter_components_by_span(dm, lo=150)
    if curve == "neg_dX_Jc":
        return C.sample_curve_column(m_solid, xa, ya_left, at, half_win=2, gap_search=14, note_prefix=f"SM Fig.4 {label} -dXJc: ")
    if curve == "Rc":
        m_dot = _filter_components_by_span(dm, lo=3, hi=18)
        res = C.sample_curve_column(m_dot, xa, ya_left, at, half_win=3, gap_search=24, note_prefix=f"SM Fig.4 {label} Rc: ")
        for i, r in enumerate(res):
            if r.get("y") is None:
                res[i] = C.sample_curve_column(
                    m_solid, xa, ya_left, [at[i]], half_win=2, gap_search=12, note_prefix=f"SM Fig.4 {label} Rc(wing): "
                )[0]
        return res
    raise ValueError(f"Unknown SM Fig. 4 curve: {curve}")


# ==============================================================================
# 11. SM Fig. 5(a): Central Peak Zoom of LP3 (C, Na, Ni, V on X in [14.7, 16.6])
# ==============================================================================
def read_smfig5a(img: np.ndarray, at: list[float], curve: str) -> list[dict]:
    fr = C.find_frame(img, include_pink_right=True)
    bt = C.find_bottom_ticks(img, fr, pad_lo=5, pad_hi=22)
    lt = C.find_left_ticks(img, fr, pad_lo=5, pad_hi=22)
    rt = C.find_right_ticks(img, fr, pad_lo=5, pad_hi=22, use_pink=True)
    if len(bt) < 2 or len(lt) < 2 or len(rt) < 2:
        raise ValueError(f"SM Fig. 5(a) ticks missing: bt={bt}, lt={lt}, rt={rt}")
    xa = C.LinearAxis([(bt[0], 14.7), (bt[-1], 16.6)])
    ya_left = C.LinearAxis([(lt[0], 2.0), (lt[-1], 0.2)])
    ya_right = C.LinearAxis([(rt[0], 0.2), (rt[-1], -0.2)])

    if curve == "V":
        m = C.interior_mask(C.pink_mask(img), fr, inset=6)
        return C.sample_curve_column(m, xa, ya_right, at, half_win=2, gap_search=10, note_prefix="SM Fig.5a V: ")

    dm = C.interior_mask(C.dark_mask(img, 145.0), fr, inset=6)
    _mask_rel_rect(dm, fr, 0.0, 0.12, 0.0, 0.18)    # '(a)'
    _mask_rel_rect(dm, fr, 0.82, 1.0, 0.0, 0.20)    # legend box

    if curve == "C":
        m = _filter_components_by_span(dm, lo=150)
        return C.sample_curve_column(m, xa, ya_left, at, half_win=2, gap_search=18, note_prefix="SM Fig.5a C: ")
    if curve == "Na":
        m = _filter_components_by_span(dm, lo=15, hi=50)
        return C.sample_curve_column(m, xa, ya_left, at, half_win=3, gap_search=36, note_prefix="SM Fig.5a Na: ")
    if curve == "Ni":
        m = _filter_components_by_span(dm, lo=3, hi=14)
        return C.sample_curve_column(m, xa, ya_left, at, half_win=3, gap_search=28, note_prefix="SM Fig.5a Ni: ")
    raise ValueError(f"Unknown SM Fig. 5(a) curve: {curve}")


# ==============================================================================
# 12. SM Fig. 5(b), (c), (d): Advective (pink) & Diffusive (gray) Flux Divergences
# ==============================================================================
def read_smfig5_flux(img: np.ndarray, at: list[float], y_top: float, y_bot: float, curve: str, panel: str) -> list[dict]:
    fr = C.find_frame(img, include_pink_right=False)
    bt = C.find_bottom_ticks(img, fr, pad_lo=6, pad_hi=28)
    lt = C.find_left_ticks(img, fr, pad_lo=6, pad_hi=28)
    if len(bt) < 2 or len(lt) < 2:
        raise ValueError(f"SM Fig. 5({panel}) ticks missing: bt={bt}, lt={lt}")
    xa = C.LinearAxis([(bt[0], 14.7), (bt[-1], 16.6)])
    ya = C.LinearAxis([(lt[0], y_top), (lt[-1], y_bot)])

    R = img[:, :, 0].astype(np.int16)
    G = img[:, :, 1].astype(np.int16)
    B = img[:, :, 2].astype(np.int16)

    if curve == "advective":
        m = C.interior_mask((R - G > 28) & (R > 160) & (B > 60), fr, inset=6)
        _mask_rel_rect(m, fr, 0.39, 0.61, 0.655, 0.695)
        return C.sample_curve_column(m, xa, ya, at, half_win=3, gap_search=20, note_prefix=f"SM Fig.5{panel} adv: ")

    if curve == "diffusive":
        gray = (R + G + B) / 3.0
        spread = np.maximum(np.maximum(R, G), B) - np.minimum(np.minimum(R, G), B)
        m = C.interior_mask((gray >= 95.0) & (gray <= 170.0) & (spread < 18.0), fr, inset=6)
        m_faded = C.interior_mask((gray >= 170.0) & (gray <= 225.0) & (spread < 15.0), fr, inset=6)
        m_faded = ndimage.binary_opening(m_faded, structure=np.ones((7, 7), bool))
        m = m | m_faded
        _mask_rel_rect(m, fr, 0.0, 0.20, 0.0, 0.17)  # panel label (b)/(c)/(d)
        _mask_rel_rect(m, fr, 0.39, 0.61, 0.83, 0.87)
        _mask_rel_rect(m, fr, 0.62, 0.96, 0.62, 0.91)
        return C.sample_curve_column(m, xa, ya, at, half_win=3, gap_search=24, note_prefix=f"SM Fig.5{panel} diff: ")

    raise ValueError(f"Unknown SM Fig. 5({panel}) curve: {curve}")
