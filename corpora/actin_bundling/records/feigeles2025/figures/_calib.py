"""Shared axis-calibration, curve-sampling, marker-splitting, and bar/error-bar helpers for feigeles2025."""

from __future__ import annotations

import math
import numpy as np
from scipy import ndimage
from boxgraph.digitize import Axis


def clusters(indices: np.ndarray, gap: int = 3) -> list[tuple[float, int, int]]:
    """Group sorted integer indices into contiguous clusters (center, start, end)."""
    if len(indices) == 0:
        return []
    splits = np.where(np.diff(indices) > gap)[0] + 1
    out = []
    for g in np.split(indices, splits):
        out.append((float(g.mean()) + 0.5, int(g[0]), int(g[-1])))
    return out


def runs_1d(mask_1d: np.ndarray, gap: int = 1) -> list[tuple[int, int]]:
    """Return (start, end_inclusive) runs of True in a 1D boolean array."""
    idx = np.nonzero(mask_1d)[0]
    if len(idx) == 0:
        return []
    splits = np.where(np.diff(idx) > gap + 1)[0] + 1
    return [(int(g[0]), int(g[-1])) for g in np.split(idx, splits)]


def verify_axis(scale: str, marks: list[tuple[float, float]], max_residual_px: float = 2.0) -> Axis:
    """Build an Axis and verify its calibration residual is within tolerance."""
    ax = Axis(scale, marks)
    if ax.residual_px > max_residual_px:
        raise RuntimeError(
            f"Axis calibration residual {ax.residual_px:.2f} px exceeds {max_residual_px} px: {marks}"
        )
    return ax


def find_box_spines(
    image: np.ndarray,
    thresh: float = 140.0,
    row_frac: float = 0.45,
    col_frac: float = 0.45,
    row_range: tuple[int, int] | None = None,
    col_range: tuple[int, int] | None = None,
):
    """Find long horizontal and vertical spines in an image or subregion.

    Returns (dark, gray, h_lines, v_lines) where each line is (center, start, end).
    """
    gray = image.astype(float).mean(axis=2)
    dark = gray < thresh
    H, W = dark.shape
    r0, r1 = row_range if row_range is not None else (0, H)
    c0, c1 = col_range if col_range is not None else (0, W)
    sub = dark[r0:r1, c0:c1]

    min_h = row_frac * sub.shape[1]
    min_v = col_frac * sub.shape[0]
    h_cand = np.where(sub.sum(axis=1) >= min_h)[0]
    v_cand = np.where(sub.sum(axis=0) >= min_v)[0]

    h_idx = []
    for r in h_cand:
        rs = runs_1d(sub[r], gap=3)
        if rs and max(b - a + 1 for a, b in rs) >= min_h:
            h_idx.append(int(r) + r0)

    v_idx = []
    for c in v_cand:
        rs = runs_1d(sub[:, c], gap=3)
        if rs and max(b - a + 1 for a, b in rs) >= min_v:
            v_idx.append(int(c) + c0)

    return dark, gray, clusters(np.array(h_idx, dtype=int), gap=3), clusters(np.array(v_idx, dtype=int), gap=3)


def find_ticks(
    mask: np.ndarray,
    spine: tuple[float, int, int],
    span: tuple[int, int],
    orientation: str,
    side: int,
    band: tuple[int, int] = (8, 32),
    frac: float = 0.55,
    max_width: int = 25,
    gap: int = 2,
) -> list[float]:
    """Find tick positions protruding from a spine in direction `side` (-1 or +1)."""
    _, s0, s1 = spine
    a, b = span
    b0, b1 = band
    H, W = mask.shape
    if orientation == "y":
        if side < 0:
            c0, c1 = max(0, s0 - b1), max(0, s0 - b0)
        else:
            c0, c1 = min(W, s1 + 1 + b0), min(W, s1 + 1 + b1)
        prof = mask[max(0, a) : min(H, b), c0:c1].mean(axis=1)
        cs = clusters(np.where(prof >= frac)[0] + max(0, a), gap=gap)
    elif orientation == "x":
        if side > 0:
            r0, r1 = min(H, s1 + 1 + b0), min(H, s1 + 1 + b1)
        else:
            r0, r1 = max(0, s0 - b1), max(0, s0 - b0)
        prof = mask[r0:r1, max(0, a) : min(W, b)].mean(axis=0)
        cs = clusters(np.where(prof >= frac)[0] + max(0, a), gap=gap)
    else:
        raise ValueError(f"Unknown orientation {orientation!r}")
    return [c[0] for c in cs if (c[2] - c[1] + 1) <= max_width]


def match_ticks_regular(candidates: list[float], n_expected: int, tol: float = 6.0) -> list[float]:
    """Select the n_expected most evenly spaced tick positions from candidates."""
    c = sorted(candidates)
    if len(c) == n_expected:
        return c
    best = None
    best_err = float("inf")
    for i in range(len(c)):
        for j in range(i + 1, len(c)):
            step = (c[j] - c[i]) / (n_expected - 1)
            if step < 20.0:
                continue
            pts = [c[i]]
            ok = True
            for k in range(1, n_expected):
                target = c[i] + k * step
                near = min(c, key=lambda x: abs(x - target))
                if abs(near - target) > tol:
                    ok = False
                    break
                pts.append(near)
            if ok and len(set(pts)) == n_expected:
                xs = np.arange(n_expected, dtype=float)
                ps = np.array(pts, dtype=float)
                m, b = np.polyfit(xs, ps, 1)
                err = float(np.max(np.abs(ps - (m * xs + b))))
                if err < best_err:
                    best_err = err
                    best = pts
    if best is None:
        raise RuntimeError(f"Could not match {n_expected} regular ticks from {c}")
    return best


def sample_solid_curve_y(
    mask: np.ndarray,
    x_axis: Axis,
    y_axis: Axis,
    x_val: float,
    row_min: int,
    row_max: int,
    pick_fn=None,
    search_dx: int = 0,
    note_prefix: str = "",
) -> dict:
    """Sample y(x) on a solid curve mask at x_val within [row_min, row_max]."""
    c_f = float(x_axis.to_pixel(x_val))
    c = int(round(c_f))
    H, W = mask.shape
    c = max(0, min(W - 1, c))

    def get_runs(col):
        if col < 0 or col >= W:
            return []
        return [(r0 + row_min, r1 + row_min) for r0, r1 in runs_1d(mask[row_min : row_max + 1, col])]

    rs = get_runs(c)
    used_c = c
    if not rs and search_dx > 0:
        for d in range(1, search_dx + 1):
            for cc in (c - d, c + d):
                rs = get_runs(cc)
                if rs:
                    used_c = cc
                    break
            if rs:
                break
    if not rs:
        return {"x": float(x_val), "y": None, "found": False, "note": f"no curve pixel at col {c}"}

    if pick_fn is not None:
        seg = pick_fn(rs, used_c)
    else:
        seg = rs[0] if len(rs) == 1 else max(rs, key=lambda t: t[1] - t[0])
    if seg is None:
        return {"x": float(x_val), "y": None, "found": False, "note": f"no matching run at col {used_c}"}

    r0, r1 = seg
    r_mid = 0.5 * (r0 + r1) + 0.5
    thick = float(r1 - r0 + 1)

    pts = []
    for dc in (-5, -3, -1, 0, 1, 3, 5):
        cc = used_c + dc
        s_rs = get_runs(cc)
        if s_rs:
            near = min(s_rs, key=lambda t: abs(0.5 * (t[0] + t[1]) + 0.5 - r_mid))
            nm = 0.5 * (near[0] + near[1]) + 0.5
            if abs(nm - r_mid) <= max(45.0, 2.5 * thick):
                pts.append((float(cc) + 0.5, nm))
    slope = 0.0
    if len(pts) >= 2:
        xs_p = np.array([p[0] for p in pts])
        ys_p = np.array([p[1] for p in pts])
        if np.ptp(xs_p) > 0:
            slope = float(np.polyfit(xs_p, ys_p, 1)[0])

    r_est = r_mid + slope * (c_f - (used_c + 0.5))
    half_px = max(0.5 * thick, 1.5) + 0.5 + y_axis.residual_px + abs(slope) * (0.5 + x_axis.residual_px)
    y_val = float(y_axis.to_value(r_est))
    y_a = float(y_axis.to_value(r_est - half_px))
    y_b = float(y_axis.to_value(r_est + half_px))
    lo, hi = min(y_a, y_b), max(y_a, y_b)
    return {
        "x": float(x_val),
        "y": y_val,
        "lo": lo,
        "hi": hi,
        "px": [c_f, float(r_est)],
        "px_lo": float(r_est - half_px),
        "px_hi": float(r_est + half_px),
        "note": f"{note_prefix}cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px, thick={thick:.0f}px",
    }


def sample_curve_local_fit(
    mask: np.ndarray,
    x_axis: Axis,
    y_axis: Axis,
    x_val: float,
    row_min: int,
    row_max: int,
    window_px: int = 55,
    deg: int = 2,
    note_prefix: str = "",
) -> dict:
    """Sample a smooth solid curve at x_val, interpolating smoothly if a scatter marker occludes the exact column."""
    # First check direct column; only use direct if unobstructed (single unsplit vertical run)
    c_f = float(x_axis.to_pixel(x_val))
    c = int(round(c_f))
    H, W = mask.shape
    direct_rs = runs_1d(mask[row_min : row_max + 1, c], gap=2) if 0 <= c < W else []
    if len(direct_rs) == 1:
        direct = sample_solid_curve_y(mask, x_axis, y_axis, x_val, row_min, row_max, search_dx=0, note_prefix=note_prefix)
        if direct.get("y") is not None:
            return direct

    c_lo = max(0, c - window_px)
    c_hi = min(W - 1, c + window_px)
    xs_l, ys_l, ths = [], [], []
    for cc in range(c_lo, c_hi + 1):
        rs = runs_1d(mask[row_min : row_max + 1, cc], gap=2)
        if len(rs) != 1:
            continue
        seg = rs[0]
        xs_l.append(float(cc) + 0.5)
        ys_l.append(0.5 * (seg[0] + seg[1]) + 0.5 + row_min)
        ths.append(float(seg[1] - seg[0] + 1))
    if len(xs_l) < 4:
        return {"x": float(x_val), "y": None, "found": False, "note": f"no curve pixels within +/-{window_px}px of col {c}"}

    xs = np.array(xs_l) - c_f
    ys = np.array(ys_l)
    w = np.exp(-0.5 * (xs / (0.5 * window_px)) ** 2)
    eff_deg = min(deg, 1 if (np.ptp(xs) < 20 or len(xs) < 8) else deg)
    coeffs = np.polyfit(xs, ys, eff_deg, w=np.sqrt(w))
    r_est = float(np.polyval(coeffs, 0.0))
    slope = float(np.polyval(np.polyder(coeffs), 0.0))
    fit_rms = float(np.sqrt(np.average((ys - np.polyval(coeffs, xs)) ** 2, weights=w)))
    thick = float(np.median(ths))
    dist_near = float(np.min(np.abs(xs)))
    half_px = (
        max(0.5 * thick, 1.5)
        + 0.5
        + fit_rms
        + min(4.0, 0.05 * dist_near)
        + y_axis.residual_px
        + abs(slope) * (0.5 + x_axis.residual_px)
    )
    y_val = float(y_axis.to_value(r_est))
    ya_, yb_ = float(y_axis.to_value(r_est - half_px)), float(y_axis.to_value(r_est + half_px))
    return {
        "x": float(x_val),
        "y": y_val,
        "lo": min(ya_, yb_),
        "hi": max(ya_, yb_),
        "px": [c_f, r_est],
        "px_lo": float(r_est - half_px),
        "px_hi": float(r_est + half_px),
        "note": f"{note_prefix}bridged occluded col (d={dist_near:.0f}px, rms={fit_rms:.2f}px), cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
    }


def extract_disk_centers(
    mask: np.ndarray,
    radius_px: float,
    min_area_frac: float = 0.35,
    min_sep_frac: float = 0.65,
    peak_frac: float = 0.45,
) -> list[tuple[float, float, int]]:
    """Find centers (col, row, support_area) of filled circular markers, splitting touching pairs/clusters.

    Uses connected components for isolated disks and local-maxima on smoothed Euclidean
    distance transform (plus watershed assignment) for merged clusters.
    """
    single_area = math.pi * radius_px * radius_px
    min_area = max(12, int(round(min_area_frac * single_area)))
    min_sep = max(3, int(round(min_sep_frac * radius_px)))

    lab, _ = ndimage.label(mask)
    centers: list[tuple[float, float, int]] = []
    for i, sl in enumerate(ndimage.find_objects(lab), 1):
        if sl is None:
            continue
        sub = lab[sl] == i
        area = int(sub.sum())
        if area < min_area:
            continue
        r0, c0 = sl[0].start, sl[1].start
        h, w = sub.shape
        # Check if isolated single disk: compact bounding box and area ~ single_area
        if area <= 1.45 * single_area and max(h, w) <= 2.65 * radius_px:
            rr, cc = np.nonzero(sub)
            centers.append((float(cc.mean()) + c0 + 0.5, float(rr.mean()) + r0 + 0.5, area))
            continue

        # Merged cluster: pad and run smoothed distance transform + local-maxima peak detection
        pad = int(math.ceil(radius_px)) + 4
        padded = np.pad(sub, pad, mode="constant", constant_values=False)
        dist = ndimage.distance_transform_edt(padded)
        dist_s = ndimage.gaussian_filter(dist, sigma=max(1.0, 0.12 * radius_px))
        # Find local maxima in footprint of radius min_sep
        yy, xx = np.ogrid[-min_sep : min_sep + 1, -min_sep : min_sep + 1]
        fp = (xx * xx + yy * yy) <= min_sep * min_sep
        is_max = (dist_s == ndimage.maximum_filter(dist_s, footprint=fp)) & (dist >= peak_frac * radius_px)
        pr, pc = np.nonzero(is_max)
        # Sort peaks by distance descending and greedily enforce min_sep
        order = np.argsort(-dist_s[pr, pc])
        kept: list[tuple[float, float]] = []
        for idx in order:
            r_p, c_p = float(pr[idx]), float(pc[idx])
            if all(math.hypot(r_p - kr, c_p - kc) >= min_sep for kr, kc in kept):
                kept.append((r_p, c_p))
        if not kept:
            rr, cc = np.nonzero(sub)
            centers.append((float(cc.mean()) + c0 + 0.5, float(rr.mean()) + r0 + 0.5, area))
            continue
        # Refine each peak by centroid of Voronoi-assigned disk pixels within 1.15*radius_px
        rr_all, cc_all = np.nonzero(padded)
        dmat = np.stack([np.hypot(rr_all - kr, cc_all - kc) for kr, kc in kept], axis=0)
        assign = np.argmin(dmat, axis=0)
        for k, (kr, kc) in enumerate(kept):
            sel = (assign == k) & (dmat[k] <= 1.25 * radius_px)
            if sel.sum() >= max(8, int(0.2 * single_area)):
                # Weighted blend of EDT peak (unbiased for partial overlap) and core centroid
                core = (assign == k) & (dmat[k] <= 0.75 * radius_px)
                if core.sum() >= 6:
                    r_ref = float(rr_all[core].mean()) - pad + r0 + 0.5
                    c_ref = float(cc_all[core].mean()) - pad + c0 + 0.5
                else:
                    r_ref = kr - pad + r0 + 0.5
                    c_ref = kc - pad + c0 + 0.5
                centers.append((c_ref, r_ref, int(sel.sum())))
    return centers


def split_disks_xor(
    mask: np.ndarray,
    radius_px: float,
    min_new_frac: float = 0.12,
    min_sep_frac: float = 0.18,
    min_xor_gain_frac: float = 0.05,
    max_per_cluster: int = 16,
) -> list[tuple[float, float, int]]:
    """Decompose isolated and overlapping filled circular disks via greedy FFT matching pursuit + XOR refinement."""
    from scipy.signal import fftconvolve

    single_area = math.pi * radius_px * radius_px
    min_area = max(12, int(round(0.32 * single_area)))
    R_int = int(math.ceil(radius_px))
    yy, xx = np.mgrid[-R_int : R_int + 1, -R_int : R_int + 1]
    disk_tpl = ((xx * xx + yy * yy) <= radius_px * radius_px).astype(np.float32)

    def paint_union(shape: tuple[int, int], pts: list[list[float]], r: float) -> np.ndarray:
        H, W = shape
        canvas = np.zeros((H, W), bool)
        Y, X = np.ogrid[:H, :W]
        r2 = r * r
        for cr, cc in pts:
            canvas |= ((Y - cr) ** 2 + (X - cc) ** 2) <= r2
        return canvas

    def refine(padded: np.ndarray, placed: list[list[float]], r: float) -> int:
        H, W = padded.shape
        Y, X = np.ogrid[:H, :W]
        r2 = r * r
        R_box = int(math.ceil(r)) + 4
        # Count buffer: how many placed disks cover each pixel
        cov_cnt = np.zeros((H, W), dtype=np.int16)
        for cr, cc in placed:
            r0 = max(0, int(cr) - R_box)
            r1 = min(H, int(cr) + R_box + 1)
            c0 = max(0, int(cc) - R_box)
            c1 = min(W, int(cc) + R_box + 1)
            dy = np.arange(r0, r1)[:, None] - cr
            dx = np.arange(c0, c1)[None, :] - cc
            cov_cnt[r0:r1, c0:c1] += (dy * dy + dx * dx) <= r2

        def disk_patch(cr: float, cc: float, r0: int, r1: int, c0: int, c1: int) -> np.ndarray:
            dy = np.arange(r0, r1)[:, None] - cr
            dx = np.arange(c0, c1)[None, :] - cc
            return (dy * dy + dx * dx) <= r2

        for step in (2.0, 1.0, 0.5):
            for _it in range(10):
                improved = False
                for k in range(len(placed)):
                    cr, cc = placed[k]
                    r0 = max(0, int(cr) - R_box - int(math.ceil(step)))
                    r1 = min(H, int(cr) + R_box + int(math.ceil(step)) + 1)
                    c0 = max(0, int(cc) - R_box - int(math.ceil(step)))
                    c1 = min(W, int(cc) + R_box + int(math.ceil(step)) + 1)
                    cur_m = disk_patch(cr, cc, r0, r1, c0, c1)
                    other_cov = (cov_cnt[r0:r1, c0:c1] - cur_m) > 0
                    sub_p = padded[r0:r1, c0:c1]
                    best_local = int(((other_cov | cur_m) ^ sub_p).sum())
                    best_pos = (cr, cc)
                    best_m = cur_m
                    for dr, dc in (
                        (step, 0), (-step, 0), (0, step), (0, -step),
                        (step, step), (step, -step), (-step, step), (-step, -step),
                    ):
                        ncr, ncc = cr + dr, cc + dc
                        nm = disk_patch(ncr, ncc, r0, r1, c0, c1)
                        err_l = int(((other_cov | nm) ^ sub_p).sum())
                        if err_l < best_local:
                            best_local = err_l
                            best_pos = (ncr, ncc)
                            best_m = nm
                    if best_pos != (cr, cc):
                        cov_cnt[r0:r1, c0:c1] += best_m.astype(np.int16) - cur_m.astype(np.int16)
                        placed[k] = [best_pos[0], best_pos[1]]
                        improved = True
                if not improved:
                    break
        return int(((cov_cnt > 0) ^ padded).sum())

    lab, _ = ndimage.label(mask)
    centers: list[tuple[float, float, int]] = []
    for i, sl in enumerate(ndimage.find_objects(lab), 1):
        if sl is None:
            continue
        sub = lab[sl] == i
        area = int(sub.sum())
        if area < min_area:
            continue
        r0, c0 = sl[0].start, sl[1].start
        h, w = sub.shape
        if area <= 1.25 * single_area and max(h, w) <= 2.30 * radius_px:
            rr, cc = np.nonzero(sub)
            centers.append((float(cc.mean()) + c0 + 0.5, float(rr.mean()) + r0 + 0.5, 1))
            continue

        pad = R_int + 4
        padded = np.pad(sub, pad, mode="constant", constant_values=False)
        H, W = padded.shape
        Y, X = np.ogrid[:H, :W]
        outside = (~padded).astype(np.float32)
        spill_map = fftconvolve(outside, disk_tpl, mode="same")
        placed: list[list[float]] = []
        covered = np.zeros_like(padded)
        min_sep = min_sep_frac * radius_px
        cur_err = int(padded.sum())

        for _ in range(max_per_cluster):
            uncov = (padded & ~covered).astype(np.float32)
            if uncov.sum() < 0.08 * single_area:
                break
            gain = fftconvolve(uncov, disk_tpl, mode="same")
            score = gain - 1.5 * spill_map
            for pr, pc in placed:
                score[((Y - pr) ** 2 + (X - pc) ** 2) < min_sep**2] = -1e9
            br, bc = divmod(int(np.argmax(score)), W)
            if gain[br, bc] < min_new_frac * single_area or spill_map[br, bc] > 0.30 * single_area:
                break
            cand = [list(p) for p in placed] + [[float(br), float(bc)]]
            new_err = refine(padded, cand, radius_px)
            if cur_err - new_err < min_xor_gain_frac * single_area:
                break
            placed, cur_err = cand, new_err
            covered = paint_union((H, W), placed, radius_px)

        for pr, pc in placed:
            centers.append((pc - pad + c0 + 0.5, pr - pad + r0 + 0.5, len(placed)))
    return centers


def read_hist_bars_red_err(
    image: np.ndarray,
    xa: Axis,
    ya: Axis,
    ht: tuple[float, int, int],
    hb: tuple[float, int, int],
    vl0: tuple[float, int, int],
    vr0: tuple[float, int, int],
    min_bar_w: int = 35,
    note_prefix: str = "",
) -> list[dict]:
    """Read solid black histogram bars with red upper error bars (lower end hidden inside black bar -> mirrored)."""
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    r_top_lim = int(ht[2]) + 10
    r_bot_lim = int(hb[1]) - 4
    c_left = int(vl0[2]) + 25
    c_right = int(vr0[1]) - 25

    black = (R < 70) & (G < 70) & (B < 70)
    red = (R > 120) & (R - G > 50) & (R - B > 50)

    # Identify columns belonging to solid black bars standing on the bottom spine
    base_band = black[r_bot_lim - 40 : r_bot_lim, :]
    col_on = np.zeros(image.shape[1], bool)
    col_on[c_left:c_right] = base_band[:, c_left:c_right].mean(axis=0) > 0.8
    bar_segs = [(a, b) for a, b in runs_1d(col_on, gap=2) if (b - a + 1) >= min_bar_w]

    out = []
    for c0, c1 in bar_segs:
        c_mid = 0.5 * (c0 + c1) + 0.5
        # Find top edge of black bar in interior columns (trim 5 px from edges)
        ic0, ic1 = c0 + 5, c1 - 5
        tops = []
        for c in range(ic0, ic1 + 1):
            # Walk upward from r_bot_lim while black
            col_b = black[r_top_lim : r_bot_lim + 1, c]
            rs = runs_1d(col_b, gap=1)
            if rs:
                tops.append(r_top_lim + rs[-1][0])
        top_row = float(np.median(tops))
        y_val = float(ya.to_value(top_row))
        x_val = float(xa.to_value(c_mid))

        # Trace red error bar upward from top_row around c_mid +/- 25 px
        stem_c0, stem_c1 = int(round(c_mid)) - 30, int(round(c_mid)) + 31
        red_rows = np.nonzero(red[r_top_lim : int(round(top_row)) + 2, stem_c0:stem_c1].any(axis=1))[0]
        if len(red_rows) == 0:
            raise RuntimeError(f"No red error bar found above bar at x={x_val:.2f}")
        cap_r_min = r_top_lim + int(red_rows.min())
        # Refine horizontal cap center within 12 px below cap_r_min
        best_w = 0
        cap_widths = []
        for rr in range(cap_r_min, min(cap_r_min + 14, int(round(top_row)))):
            w = int(red[rr, stem_c0:stem_c1].sum())
            cap_widths.append((rr, w))
            if w > best_w:
                best_w = w
        cap_rows = [rr for rr, w in cap_widths if w >= 0.7 * best_w]
        row_hi = float(np.mean(cap_rows)) + 0.5

        half_px = 1.5 + ya.residual_px
        ya_, yb_ = float(ya.to_value(top_row - half_px)), float(ya.to_value(top_row + half_px))
        end_unc = 2.0 + ya.residual_px
        bar_hi = float(ya.to_value(row_hi))
        bar_hi_outer = float(ya.to_value(row_hi - end_unc))
        # Lower error bar hidden inside the solid black bar: mirror upper error bar
        bar_lo = 2.0 * y_val - bar_hi
        bar_lo_outer = 2.0 * y_val - bar_hi_outer
        row_lo = float(ya.to_pixel(bar_lo))

        out.append(
            {
                "x": x_val,
                "y": y_val,
                "lo": min(ya_, yb_),
                "hi": max(ya_, yb_),
                "bar_lo": bar_lo,
                "bar_hi": bar_hi,
                "bar_lo_outer": bar_lo_outer,
                "bar_hi_outer": bar_hi_outer,
                "px": [c_mid, top_row],
                "px_lo": row_hi,
                "px_hi": row_lo,
                "note": (
                    f"{note_prefix}bar w={c1-c0+1}px, upper red bar={bar_hi:.3f} (+{bar_hi-y_val:.3f}); "
                    f"lower end hidden inside black bar and mirrored; "
                    f"cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px"
                ),
            }
        )
    return out


def make_marker_reading(
    c: float,
    r: float,
    x_axis: Axis,
    y_axis: Axis,
    half_px: float,
    note: str,
) -> dict:
    """Convert pixel center (c, r) to a calibrated BoxGraph reading dict."""
    x_val = float(x_axis.to_value(c))
    y_val = float(y_axis.to_value(r))
    ya = float(y_axis.to_value(r - half_px))
    yb = float(y_axis.to_value(r + half_px))
    return {
        "x": x_val,
        "y": y_val,
        "lo": min(ya, yb),
        "hi": max(ya, yb),
        "px": [float(c), float(r)],
        "px_lo": float(r - half_px),
        "px_hi": float(r + half_px),
        "note": note,
    }
