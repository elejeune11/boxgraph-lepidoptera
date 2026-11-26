"""Calibrated image reader for S4 Figure (p. 24 / SI-5): bound crosslinker mole ratio alpha_bound:A
vs total crosslinker mole ratio alpha:A across 6 M:A ratios for unipolar (left) and apolar (right)
bundles (900 dpi).

Extracts all 6 series by fitting the straight connecting line segments between consecutive alpha:A
values, stitching across vertices starting from the well-separated [0.4, 0.6] interval (strictly
ordered top-to-bottom as M:A = 0.0225, 0.045, 0.09, 0.18, 0.225, 0.675), and locating each point's
symmetric +/- SD horizontal error-bar caps directly in the image.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _common as C

MAS = [0.0225, 0.045, 0.09, 0.18, 0.225, 0.675]
ALL_X = [0.01, 0.05, 0.1, 0.25, 0.4, 0.6, 0.8]


def _calibrate(image: np.ndarray):
    info = C.find_matlab_box(image, dark_thresh=165.0, min_tick_len=5, max_tick_len=55)
    sp = info["spines"]
    xt = C.merge_tick_set(sp["left"], info["x_ticks_top"] + info["x_ticks_bot"], sp["right"], 9, tol_px=20.0)
    yt = C.merge_tick_set(sp["top"], info["y_ticks_left"], sp["bottom"], 10, tol_px=20.0)
    xa = C.fit_axis(xt, [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], max_residual_px=2.5)
    ya = C.fit_axis(yt, [0.045, 0.04, 0.035, 0.03, 0.025, 0.02, 0.015, 0.01, 0.005, 0.0], max_residual_px=2.5)
    return info, xa, ya


def _mask_legend_and_ticks(dark: np.ndarray, box: tuple[int, int, int, int]) -> np.ndarray:
    c0, r0, c1, r1 = box
    clean = dark.copy()
    clean[:r0 + 24, :] = False
    clean[r1 - 16:, :] = False
    clean[:, :c0 + 20] = False
    clean[:, c1 - 20:] = False
    # Fast legend box detection on the right half of the plot box
    x_mid = c0 + (c1 - c0) // 2
    sub = dark[r0 + 30:r1 - 30, x_mid:c1 - 15]
    h_lines = C.find_spines_1d(sub, min_len=260, gap=3)
    for i in range(len(h_lines)):
        for j in range(i + 1, len(h_lines)):
            a, b = h_lines[i], h_lines[j]
            if (
                abs(a["start"] - b["start"]) <= 12
                and abs(a["end"] - b["end"]) <= 12
                and 200 <= (b["lo"] - a["hi"]) <= 550
                and 260 <= a["length"] <= 500
            ):
                lr0 = r0 + 30 + a["lo"] - 6
                lr1 = r0 + 30 + b["hi"] + 6
                lc0 = x_mid + min(a["start"], b["start"]) - 6
                lc1 = x_mid + max(a["end"], b["end"]) + 6
                clean[lr0:lr1, lc0:lc1] = False
                return clean
    return clean


def _extract_col_runs(mask: np.ndarray, col: int, max_run_h: int = 40) -> list[tuple[float, int]]:
    rs = np.flatnonzero(mask[:, col])
    if rs.size == 0:
        return []
    splits = np.flatnonzero(np.diff(rs) > 2)
    st = np.concatenate(([0], splits + 1))
    en = np.concatenate((splits, [rs.size - 1]))
    out = []
    for a, b in zip(st, en):
        h = int(rs[b] - rs[a] + 1)
        if 2 <= h <= max_run_h:
            out.append((0.5 * (rs[a] + rs[b]) + 0.5, h))
    return out


def _fit_segments_in_interval(clean: np.ndarray, col_L: float, col_R: float,
                              pad_L: int = 34, pad_R: int = 34,
                              max_run_h: int = 40) -> list[tuple[float, float, float, float]]:
    """Track straight line segments across [col_L + pad_L, col_R - pad_R]."""
    c_start = int(round(col_L)) + pad_L
    c_end = int(round(col_R)) - pad_R
    if c_end - c_start < 12:
        c_start = int(round(col_L)) + max(14, int(0.18 * (col_R - col_L)))
        c_end = int(round(col_R)) - max(14, int(0.18 * (col_R - col_L)))
    step = max(2, (c_end - c_start) // 28)
    cols = list(range(c_start, c_end + 1, step))

    tracks: list[list[tuple[int, float]]] = []
    for c in cols:
        runs = _extract_col_runs(clean, c, max_run_h=max_run_h)
        used = set()
        preds = []
        for ti, tr in enumerate(tracks):
            if c - tr[-1][0] > 4 * step:
                continue
            if len(tr) >= 2:
                k0 = max(0, len(tr) - 6)
                xs = np.array([p[0] for p in tr[k0:]], dtype=float)
                ys = np.array([p[1] for p in tr[k0:]], dtype=float)
                m, b = np.polyfit(xs, ys, 1)
                pred = m * c + b
            else:
                pred = tr[-1][1]
            preds.append((ti, pred))
        pairs = []
        for ti, pred in preds:
            for ri, (r_val, _) in enumerate(runs):
                err = abs(r_val - pred)
                tol = 18.0 if len(tracks[ti]) == 1 else 10.0
                if err <= tol:
                    pairs.append((err, ti, ri))
        pairs.sort()
        matched_t = set()
        for err, ti, ri in pairs:
            if ti in matched_t or ri in used:
                continue
            tracks[ti].append((c, runs[ri][0]))
            matched_t.add(ti)
            used.add(ri)
        for ri, (r_val, _) in enumerate(runs):
            if ri not in used:
                tracks.append([(c, r_val)])

    min_pts = max(4, int(0.40 * len(cols)))
    segs = []
    for tr in tracks:
        if len(tr) < min_pts:
            continue
        xs = np.array([p[0] for p in tr], dtype=float)
        ys = np.array([p[1] for p in tr], dtype=float)
        m, b = np.polyfit(xs, ys, 1)
        res = np.abs(ys - (m * xs + b))
        inliers = res <= max(2.5, float(np.median(res)) * 2.5)
        if inliers.sum() >= max(3, min_pts - 1):
            m, b = np.polyfit(xs[inliers], ys[inliers], 1)
            res_max = float(np.max(np.abs(ys[inliers] - (m * xs[inliers] + b))))
        else:
            res_max = float(np.max(res))
        segs.append((float(m * col_L + b), float(m * col_R + b), float(m), res_max))

    segs.sort(key=lambda s: 0.5 * (s[0] + s[1]))
    dedup = []
    for s in segs:
        if dedup and abs(s[0] - dedup[-1][0]) < 8.0 and abs(s[1] - dedup[-1][1]) < 8.0:
            if s[3] < dedup[-1][3]:
                dedup[-1] = s
        else:
            dedup.append(s)
    return dedup


def _detect_all_caps_in_col(dark: np.ndarray, box: tuple[int, int, int, int],
                            cx: int, is_right_border: bool,
                            cap_arm: tuple[int, int], y_ticks: list[float]) -> list[float]:
    """Find all horizontal error-bar cap center rows in column cx."""
    c0, r0, c1, r1 = box
    H, W = dark.shape
    a0, a1 = cap_arm
    r_lo = r0 + 4
    r_hi = r1 - 2

    L_arm = dark[r_lo:r_hi, max(c0, cx - a1):max(c0 + 1, cx - a0)].mean(axis=1)
    if is_right_border:
        R_arm = dark[r_lo:r_hi, min(W - 1, cx + a0):min(W, cx + a1)].mean(axis=1)
        # Error bar caps at x=0.8 extend OUTSIDE the right spine (R_arm) AND inside (L_arm),
        # whereas right-spine y-ticks extend ONLY inside (L_arm, R_arm == 0)!
        cap_score = np.minimum(L_arm, R_arm)
    else:
        R_arm = dark[r_lo:r_hi, min(c1 - 1, cx + a0):min(c1, cx + a1)].mean(axis=1)
        cap_score = np.minimum(L_arm, R_arm)

    rows = np.flatnonzero(cap_score >= 0.55) + r_lo
    caps: list[float] = []
    if rows.size:
        splits = np.flatnonzero(np.diff(rows) > 2)
        st = np.concatenate(([0], splits + 1))
        en = np.concatenate((splits, [rows.size - 1]))
        for a, b in zip(st, en):
            h = int(rows[b] - rows[a] + 1)
            # Standard cap is 4..9 px tall; two closely touching caps can merge to 12..18 px
            if 2 <= h <= 13:
                caps.append(0.5 * (rows[a] + rows[b]) + 0.5)
            elif 14 <= h <= 20:
                # Two merged caps: record top edge + 3.5 and bottom edge - 3.5
                caps.append(float(rows[a]) + 3.5)
                caps.append(float(rows[b]) - 3.5)
    return caps


def _assign_caps_at_column(dark: np.ndarray, box: tuple[int, int, int, int],
                           cx: int, col_verts: dict[float, float], is_right_border: bool,
                           cap_arm: tuple[int, int], y_ticks: list[float]) -> dict[float, tuple[float | None, float | None, str]]:
    """Assign symmetric +/- SD caps to all vertices at column cx globally (pairs first, then unclaimed singles)."""
    H, W = dark.shape
    c0, r0, c1, r1 = box
    caps = _detect_all_caps_in_col(dark, box, cx, is_right_border, cap_arm, y_ticks)
    stem_col = dark[:, max(0, cx - 3):min(W, cx + 4)].any(axis=1)

    def has_stem(r1_px: float, r2_px: float) -> bool:
        a, b = int(round(min(r1_px, r2_px))), int(round(max(r1_px, r2_px)))
        if b <= a:
            return False
        return float(stem_col[a:b + 1].mean()) >= 0.85

    # Reject caps that fall right inside any marker's solid fill (within 15 px of a vertex for non-0.0225/0.675)
    # Wait: the `x` marker (0.675) has 4 arms at ~7..10 px from center, so require min_d >= 16 px for filled markers
    def min_dist_for(ma: float) -> float:
        if ma == 0.0225:
            return 10.0
        if ma == 0.675:
            return 16.0
        return 19.0

    # Pass 1: find all valid symmetric cap pairs (above, below) for each vertex
    pair_cands = []
    for ma, rv in col_verts.items():
        md = min_dist_for(ma)
        above = [i for i, rc in enumerate(caps) if rc < rv - md and (rv - rc) <= 210.0 and has_stem(rv, rc)]
        below = [i for i, rc in enumerate(caps) if rc > rv + md and (rc - rv) <= 210.0 and has_stem(rv, rc)]
        for ia in above:
            dt = rv - caps[ia]
            for ib in below:
                db = caps[ib] - rv
                asym = abs(dt - db)
                if asym <= max(6.5, 0.09 * 0.5 * (dt + db)):
                    pair_cands.append((asym + 0.015 * (dt + db), ma, ia, ib))

    pair_cands.sort()
    assigned: dict[float, tuple[float | None, float | None, str]] = {}
    used_caps: set[int] = set()
    for _, ma, ia, ib in pair_cands:
        if ma in assigned:
            continue
        # Allow cap sharing ONLY if unclaimed caps don't exist; prefer strictly unclaimed caps first
        if ia in used_caps or ib in used_caps:
            continue
        assigned[ma] = (caps[ia], caps[ib], "both caps")
        used_caps.add(ia)
        used_caps.add(ib)

    # Pass 1b: symmetric pairs that may share one touching cap with an already-assigned neighbor
    for _, ma, ia, ib in pair_cands:
        if ma in assigned:
            continue
        if not (ia in used_caps and ib in used_caps):
            assigned[ma] = (caps[ia], caps[ib], "both caps")
            used_caps.add(ia)
            used_caps.add(ib)

    # Pass 2: for remaining vertices (where one cap is hidden inside a neighboring marker or bottom spine),
    # pick the closest UNCLAIMED cap connected by a continuous stem and mirror it symmetrically!
    for ma, rv in col_verts.items():
        if ma in assigned:
            continue
        md = min_dist_for(ma)
        cands = []
        for i, rc in enumerate(caps):
            if i in used_caps:
                continue
            d = abs(rc - rv)
            if md <= d <= 185.0 and has_stem(rv, rc):
                # Prefer caps not sitting directly inside another vertex's marker body
                near_other = any(abs(rc - rv2) < 12.0 for ma2, rv2 in col_verts.items() if ma2 != ma)
                if not near_other:
                    cands.append((d, i, rc))
        if cands:
            cands.sort()
            d, i, rc = cands[0]
            used_caps.add(i)
            if rc < rv:
                assigned[ma] = (rc, rv + d, "upper cap mirrored")
            else:
                assigned[ma] = (rv - d, rc, "lower cap mirrored")
        else:
            assigned[ma] = (None, None, "no cap")

    return assigned


def extract_panel(image: np.ndarray, polarity: str) -> tuple[dict, dict, object, object]:
    """Extract all 6 M:A curves (vertices + symmetric SD error bars) for S4 unipolar or apolar."""
    info, xa, ya = _calibrate(image)
    dark = (image[:, :, 0].astype(np.uint16) + image[:, :, 1] + image[:, :, 2]) < 450
    clean = _mask_legend_and_ticks(dark, info["box"])

    x_cols = {xv: xa.to_pixel(xv) for xv in ALL_X}

    segs_46 = _fit_segments_in_interval(clean, x_cols[0.4], x_cols[0.6], pad_L=36, pad_R=36)
    segs_68 = _fit_segments_in_interval(clean, x_cols[0.6], x_cols[0.8], pad_L=36, pad_R=36)
    segs_254 = _fit_segments_in_interval(clean, x_cols[0.25], x_cols[0.4], pad_L=36, pad_R=36)
    segs_125 = _fit_segments_in_interval(clean, x_cols[0.1], x_cols[0.25], pad_L=32, pad_R=36, max_run_h=55)
    segs_051 = _fit_segments_in_interval(clean, x_cols[0.05], x_cols[0.1], pad_L=24, pad_R=24, max_run_h=65)

    if len(segs_46) != 6:
        raise ValueError(f"S4 {polarity}: expected 6 segments on [0.4, 0.6], found {len(segs_46)}")

    verts: dict[float, dict[float, float]] = {ma: {} for ma in MAS}
    for ma, seg in zip(MAS, segs_46):
        verts[ma][0.4] = seg[0]
        verts[ma][0.6] = seg[1]

    def match_right(src_x: float, dst_x: float, segs: list):
        used = set()
        pairs = []
        for ma in MAS:
            if src_x not in verts[ma]:
                continue
            r0 = verts[ma][src_x]
            for si, s in enumerate(segs):
                pairs.append((abs(s[0] - r0), ma, si))
        pairs.sort()
        matched_ma = set()
        for err, ma, si in pairs:
            if ma in matched_ma or si in used or err > 20.0:
                continue
            verts[ma][src_x] = 0.5 * (verts[ma][src_x] + segs[si][0])
            verts[ma][dst_x] = segs[si][1]
            matched_ma.add(ma)
            used.add(si)

    def match_left(src_x: float, dst_x: float, segs: list, tol: float = 22.0):
        used = set()
        pairs = []
        for ma in MAS:
            if src_x not in verts[ma]:
                continue
            r0 = verts[ma][src_x]
            for si, s in enumerate(segs):
                pairs.append((abs(s[1] - r0), ma, si))
        pairs.sort()
        matched_ma = set()
        for err, ma, si in pairs:
            if ma in matched_ma or si in used or err > tol:
                continue
            verts[ma][src_x] = 0.5 * (verts[ma][src_x] + segs[si][1])
            verts[ma][dst_x] = segs[si][0]
            matched_ma.add(ma)
            used.add(si)

    match_right(0.6, 0.8, segs_68)
    match_left(0.4, 0.25, segs_254)
    match_left(0.25, 0.1, segs_125, tol=25.0)
    match_left(0.1, 0.05, segs_051, tol=25.0)

    cap_arm = (5, 13) if polarity == "unipolar" else (4, 10)
    bars: dict[float, dict[float, tuple[float | None, float | None, str]]] = {ma: {} for ma in MAS}
    for xv in [0.05, 0.1, 0.25, 0.4, 0.6, 0.8]:
        col_verts = {ma: verts[ma][xv] for ma in MAS if xv in verts[ma]}
        cx = int(round(x_cols[xv]))
        assigned = _assign_caps_at_column(
            dark, info["box"], cx, col_verts, is_right_border=(xv == 0.8),
            cap_arm=cap_arm, y_ticks=info["y_ticks_left"],
        )
        for ma, (rt, rb, tag) in assigned.items():
            if rt is not None and rb is not None and tag == "both caps":
                mid = 0.5 * (rt + rb)
                if abs(mid - verts[ma][xv]) <= 5.0:
                    verts[ma][xv] = 0.5 * (verts[ma][xv] + mid)
            bars[ma][xv] = (rt, rb, tag)

    return verts, bars, xa, ya


def read_series(image: np.ndarray, polarity: str, ma: float, at: list[float]) -> list[dict]:
    verts, bars, xa, ya = extract_panel(image, polarity)
    xs = at if at else sorted(verts[ma].keys())
    out = []
    for xv in xs:
        match_x = [k for k in verts[ma] if abs(k - xv) < 1e-4]
        if not match_x:
            out.append({"x": float(xv), "y": None, "note": f"S4 {polarity} M:A={ma}: x={xv} not resolved"})
            continue
        xk = match_x[0]
        rv = verts[ma][xk]
        rt, rb, tag = bars[ma][xk]
        pc = xa.to_pixel(xk)
        note = f"S4 {polarity} M:A={ma:g} at alpha:A={xk:g} ({tag}; rx={xa.residual_px:.2f}px, ry={ya.residual_px:.2f}px)"
        out.append(
            C.make_reading(
                xk, rv, ya, pc, half_px=2.0, x_axis=xa, note=note,
                bar_top_px=rt, bar_bot_px=rb, bar_err_px=2.5, clamp_y=(0.0, 0.045),
            )
        )
    return out
