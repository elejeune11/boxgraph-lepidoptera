"""Shared fast 2D affine calibration and marker/curve extraction for tang1996 figures.

Every axis calibration and data reading is derived strictly from the rendered
image pixels and the printed tick labels visible on each figure panel.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
import numpy as np
from scipy import ndimage, signal


@dataclass
class PanelCal:
    """2D affine calibration mapping image (col, row) <-> data (x, y)."""

    x_scale: str
    y_scale: str
    ax: np.ndarray  # u_x = ax[0] + ax[1]*c + ax[2]*r
    ay: np.ndarray  # u_y = ay[0] + ay[1]*r + ay[2]*c
    resid_x_px: float
    resid_y_px: float
    box: tuple[float, float, float, float]  # (c_left, r_top, c_right, r_bot)
    inner_box: tuple[int, int, int, int]    # (c0, r0, c1, r1) inside frame stroke

    def px_to_data(self, c: float, r: float) -> tuple[float, float]:
        ux = float(self.ax[0] + self.ax[1] * c + self.ax[2] * r)
        uy = float(self.ay[0] + self.ay[1] * r + self.ay[2] * c)
        x = 10.0**ux if self.x_scale == "log10" else ux
        y = 10.0**uy if self.y_scale == "log10" else uy
        return x, y

    def data_to_px(self, x: float, y: float) -> tuple[float, float]:
        ux = math.log10(x) if self.x_scale == "log10" else float(x)
        uy = math.log10(y) if self.y_scale == "log10" else float(y)
        det = self.ax[1] * self.ay[1] - self.ax[2] * self.ay[2]
        rhs_x = ux - self.ax[0]
        rhs_y = uy - self.ay[0]
        c = float((self.ay[1] * rhs_x - self.ax[2] * rhs_y) / det)
        r = float((-self.ay[2] * rhs_x + self.ax[1] * rhs_y) / det)
        return c, r

    def col_for_x(self, x: float, r: float) -> float:
        ux = math.log10(x) if self.x_scale == "log10" else float(x)
        return float((ux - self.ax[0] - self.ax[2] * r) / self.ax[1])

    def row_for_y(self, y: float, c: float) -> float:
        uy = math.log10(y) if self.y_scale == "log10" else float(y)
        return float((uy - self.ay[0] - self.ay[2] * c) / self.ay[1])

    def y_interval(self, c: float, r: float, half_px: float) -> tuple[float, float, float, float, float]:
        """Return (y, lo, hi, px_lo, px_hi) for pixel (c, r) with vertical half-error half_px."""
        _, y = self.px_to_data(c, r)
        _, y1 = self.px_to_data(c, r - half_px)
        _, y2 = self.px_to_data(c, r + half_px)
        return y, min(y1, y2), max(y1, y2), r - half_px, r + half_px

    def data_box_to_px(self, x0: float, y0: float, x1: float, y1: float) -> tuple[int, int, int, int]:
        """Convert data rectangle (x_min, y_min, x_max, y_max) to pixel bbox (c0, r0, c1, r1)."""
        ca, ra = self.data_to_px(x0, y1)
        cb, rb = self.data_to_px(x1, y0)
        return (int(min(ca, cb)), int(min(ra, rb)), int(max(ca, cb)), int(max(ra, rb)))


def make_reading(cal: PanelCal, c: float, r: float, half_px: float, note: str = "") -> dict:
    x, _ = cal.px_to_data(c, r)
    if abs(x) < 1e-10:
        x = 0.0
    y, lo, hi, px_lo, px_hi = cal.y_interval(c, r, half_px)
    return {
        "x": x,
        "y": y,
        "lo": lo,
        "hi": hi,
        "px": [float(c), float(r)],
        "px_lo": float(px_lo),
        "px_hi": float(px_hi),
        "note": note,
    }


def _to_dark(image: np.ndarray, thr: int = 160) -> np.ndarray:
    """Fast grayscale dark-ink mask (uint8, zero float64 allocation)."""
    if image.ndim == 3:
        return image[:, :, 0] < thr
    return image < thr


def _find_frame(dark: np.ndarray):
    H, W = dark.shape
    rs = dark.sum(axis=1)
    long_r = np.where(rs > 0.45 * W)[0]
    groups_r = np.split(long_r, np.where(np.diff(long_r) > 8)[0] + 1)
    top_g, bot_g = groups_r[0], groups_r[-1]
    r0, r1 = float(top_g.mean()), float(bot_g.mean())

    cs = dark[int(r0):int(r1) + 1, :].sum(axis=0)
    long_c = np.where(cs > 0.65 * (r1 - r0))[0]
    groups_c = np.split(long_c, np.where(np.diff(long_c) > 8)[0] + 1)
    left_g, right_g = groups_c[0], groups_c[-1]
    c0, c1 = float(left_g.mean()), float(right_g.mean())

    c_in0, r_in0, c_in1, r_in1 = int(left_g[-1]), int(top_g[-1]), int(right_g[0]), int(bot_g[0])

    def fit_edge_h(r_est):
        y_lo, y_hi = max(0, int(r_est) - 16), min(H, int(r_est) + 17)
        band = dark[y_lo:y_hi, :]
        cols = np.arange(c_in0 + 40, c_in1 - 40, 20)
        pts_c, pts_r = [], []
        for c in cols:
            idx = np.nonzero(band[:, c])[0]
            if 2 <= len(idx) <= 18:
                pts_c.append(c)
                pts_r.append(y_lo + float(idx.mean()))
        pts_c, pts_r = np.array(pts_c, float), np.array(pts_r, float)
        m, b = np.polyfit(pts_c, pts_r, 1)
        keep = np.abs(pts_r - (m * pts_c + b)) < 3.0
        if keep.sum() >= 5:
            m, b = np.polyfit(pts_c[keep], pts_r[keep], 1)
        return float(m), float(b)

    def fit_edge_v(c_est):
        x_lo, x_hi = max(0, int(c_est) - 16), min(W, int(c_est) + 17)
        band = dark[:, x_lo:x_hi]
        rows = np.arange(r_in0 + 40, r_in1 - 40, 20)
        pts_r, pts_c = [], []
        for r in rows:
            idx = np.nonzero(band[r, :])[0]
            if 2 <= len(idx) <= 18:
                pts_r.append(r)
                pts_c.append(x_lo + float(idx.mean()))
        pts_r, pts_c = np.array(pts_r, float), np.array(pts_c, float)
        m, b = np.polyfit(pts_r, pts_c, 1)
        keep = np.abs(pts_c - (m * pts_r + b)) < 3.0
        if keep.sum() >= 5:
            m, b = np.polyfit(pts_r[keep], pts_c[keep], 1)
        return float(m), float(b)

    e_top = fit_edge_h(r0)
    e_bot = fit_edge_h(r1)
    e_left = fit_edge_v(c0)
    e_right = fit_edge_v(c1)
    return (c0, r0, c1, r1), (c_in0, r_in0, c_in1, r_in1), (e_top, e_bot, e_left, e_right)


def _scan_ticks_h(dark: np.ndarray, m_line: float, b_line: float, direction: int, c_lo: int, c_hi: int):
    cols = np.arange(c_lo, c_hi)
    r_base = np.rint(m_line * cols + b_line).astype(int) + direction * 6
    offs = np.arange(95)[:, None] * direction
    rows = np.clip(r_base[None, :] + offs, 0, dark.shape[0] - 1)
    strip = dark[rows, cols[None, :]]
    starts = strip[:4, :].any(axis=0)
    bridged = strip.copy()
    bridged[:-1, :] |= strip[1:, :]
    first_false = np.argmax(~bridged, axis=0)
    all_true = bridged.all(axis=0)
    depths = np.where(starts, np.where(all_true, 95, first_false), 0)

    idx = np.nonzero((depths >= 16) & (depths <= 88))[0]
    if len(idx) == 0:
        return []
    out = []
    for gr in np.split(idx, np.where(np.diff(idx) > 4)[0] + 1):
        if 2 <= len(gr) <= 18:
            c_pos = float(gr.mean()) + c_lo
            r_pos = m_line * c_pos + b_line
            out.append((c_pos, r_pos, int(depths[gr].max())))
    return out


def _scan_ticks_v(dark: np.ndarray, m_line: float, b_line: float, direction: int, r_lo: int, r_hi: int):
    rows = np.arange(r_lo, r_hi)
    c_base = np.rint(m_line * rows + b_line).astype(int) + direction * 6
    offs = np.arange(95)[:, None] * direction
    cols = np.clip(c_base[None, :] + offs, 0, dark.shape[1] - 1)
    strip = dark[rows[None, :], cols]
    starts = strip[:4, :].any(axis=0)
    bridged = strip.copy()
    bridged[:-1, :] |= strip[1:, :]
    first_false = np.argmax(~bridged, axis=0)
    all_true = bridged.all(axis=0)
    depths = np.where(starts, np.where(all_true, 95, first_false), 0)

    idx = np.nonzero((depths >= 16) & (depths <= 88))[0]
    if len(idx) == 0:
        return []
    out = []
    for gr in np.split(idx, np.where(np.diff(idx) > 4)[0] + 1):
        if 2 <= len(gr) <= 18:
            r_pos = float(gr.mean()) + r_lo
            c_pos = m_line * r_pos + b_line
            out.append((c_pos, r_pos, int(depths[gr].max())))
    return out


def _u(v: float, scale: str) -> float:
    return math.log10(v) if scale == "log10" else float(v)


def _match_regular_ticks(cands: list[tuple[float, float, int]], vals: list[float], scale: str, major_min_depth: int = 48):
    maj = [t for t in cands if t[2] >= major_min_depth]
    n = len(vals)
    if len(maj) < n:
        raise ValueError(f"Expected {n} major ticks, found {len(maj)}: {maj}")
    us = np.array([_u(v, scale) for v in vals], float)
    if len(maj) == n:
        return list(zip(maj, vals))
    from itertools import combinations
    best_err = 1e9
    best_sub = None
    for sub in combinations(maj, n):
        cs = np.array([t[0] for t in sub], float)
        rs = np.array([t[1] for t in sub], float)
        ps = cs if (cs.max() - cs.min()) > (rs.max() - rs.min()) else rs
        b, a = np.polyfit(ps, us, 1)
        err = float(np.max(np.abs(ps - (us - a) / b)))
        if err < best_err:
            best_err = err
            best_sub = sub
    return list(zip(best_sub, vals))


def calibrate_panel(
    image: np.ndarray,
    *,
    x_scale: str,
    y_scale: str,
    x_major: list[float],
    y_major: list[float],
    x_all: list[float],
    y_all: list[float],
    x_edges: tuple[float | None, float | None] = (None, None),
    y_edges: tuple[float | None, float | None] = (None, None),
    x_sides: tuple[str, ...] = ("bot", "top"),
    y_sides: tuple[str, ...] = ("left", "right"),
    x_match_by_all: bool = False,
) -> PanelCal:
    """Calibrate a framed plot panel in 2D strictly from detected frame edges and inward ticks."""
    dark = _to_dark(image)
    box, inner, (e_top, e_bot, e_left, e_right) = _find_frame(dark)
    c_in0, r_in0, c_in1, r_in1 = inner

    bt = _scan_ticks_h(dark, e_bot[0], e_bot[1], -1, c_in0 + 15, c_in1 - 15)
    tt = _scan_ticks_h(dark, e_top[0], e_top[1], +1, c_in0 + 15, c_in1 - 15)
    lt = _scan_ticks_v(dark, e_left[0], e_left[1], +1, r_in0 + 15, r_in1 - 15)
    rt = _scan_ticks_v(dark, e_right[0], e_right[1], -1, r_in0 + 15, r_in1 - 15)

    def intersect_hv(e_h, e_v):
        m_h, b_h = e_h
        m_v, b_v = e_v
        c = (m_v * b_h + b_v) / (1.0 - m_v * m_h)
        r = m_h * c + b_h
        return float(c), float(r)

    tl = intersect_hv(e_top, e_left)
    tr = intersect_hv(e_top, e_right)
    bl = intersect_hv(e_bot, e_left)
    br = intersect_hv(e_bot, e_right)

    x_pts = []
    side_map_x = {"bot": (bt, bl, br), "top": (tt, tl, tr)}
    for side in x_sides:
        cands, p_left, p_right = side_map_x[side]
        if x_match_by_all:
            if len(cands) != len(x_all):
                raise ValueError(f"Expected {len(x_all)} x-ticks on {side}, got {len(cands)}: {cands}")
            for t, v in zip(cands, x_all):
                x_pts.append((t[0], t[1], _u(v, x_scale)))
            if x_edges[0] is not None:
                x_pts.append((p_left[0], p_left[1], _u(x_edges[0], x_scale)))
            if x_edges[1] is not None:
                x_pts.append((p_right[0], p_right[1], _u(x_edges[1], x_scale)))
            continue
        seed = [(t[0], t[1], _u(v, x_scale)) for t, v in _match_regular_ticks(cands, x_major, x_scale)]
        if x_edges[0] is not None:
            seed.append((p_left[0], p_left[1], _u(x_edges[0], x_scale)))
        if x_edges[1] is not None:
            seed.append((p_right[0], p_right[1], _u(x_edges[1], x_scale)))
        sc = np.array([p[0] for p in seed], float)
        su = np.array([p[2] for p in seed], float)
        b1, a1 = np.polyfit(sc, su, 1)
        for v in x_all:
            uv = _u(v, x_scale)
            c_pred = (uv - a1) / b1
            if cands:
                best = min(cands, key=lambda t: abs(t[0] - c_pred))
                if abs(best[0] - c_pred) <= 6.0:
                    x_pts.append((best[0], best[1], uv))
        if x_edges[0] is not None:
            x_pts.append((p_left[0], p_left[1], _u(x_edges[0], x_scale)))
        if x_edges[1] is not None:
            x_pts.append((p_right[0], p_right[1], _u(x_edges[1], x_scale)))

    y_pts = []
    side_map_y = {"left": (lt, tl, bl), "right": (rt, tr, br)}
    for side in y_sides:
        cands, p_top, p_bot = side_map_y[side]
        seed = [(t[0], t[1], _u(v, y_scale)) for t, v in _match_regular_ticks(cands, y_major, y_scale)]
        if y_edges[0] is not None:
            seed.append((p_top[0], p_top[1], _u(y_edges[0], y_scale)))
        if y_edges[1] is not None:
            seed.append((p_bot[0], p_bot[1], _u(y_edges[1], y_scale)))
        sr = np.array([p[1] for p in seed], float)
        su = np.array([p[2] for p in seed], float)
        b1, a1 = np.polyfit(sr, su, 1)
        for v in y_all:
            uv = _u(v, y_scale)
            r_pred = (uv - a1) / b1
            if cands:
                best = min(cands, key=lambda t: abs(t[1] - r_pred))
                if abs(best[1] - r_pred) <= 6.0:
                    y_pts.append((best[0], best[1], uv))
        if y_edges[0] is not None:
            y_pts.append((p_top[0], p_top[1], _u(y_edges[0], y_scale)))
        if y_edges[1] is not None:
            y_pts.append((p_bot[0], p_bot[1], _u(y_edges[1], y_scale)))

    X_arr = np.array(x_pts, float)
    if len(x_sides) >= 2:
        Ax = np.column_stack([np.ones(len(X_arr)), X_arr[:, 0], X_arr[:, 1]])
        ax, _, _, _ = np.linalg.lstsq(Ax, X_arr[:, 2], rcond=None)
    else:
        m_v = 0.5 * (e_left[0] + e_right[0])
        c_corr = X_arr[:, 0] - m_v * X_arr[:, 1]
        b1, a1 = np.polyfit(c_corr, X_arr[:, 2], 1)
        ax = np.array([a1, b1, -b1 * m_v], float)
    ux_pred = ax[0] + ax[1] * X_arr[:, 0] + ax[2] * X_arr[:, 1]
    resid_x_px = float(np.max(np.abs((X_arr[:, 2] - ux_pred) / ax[1])))

    Y_arr = np.array(y_pts, float)
    if len(y_sides) >= 2:
        Ay = np.column_stack([np.ones(len(Y_arr)), Y_arr[:, 1], Y_arr[:, 0]])
        ay, _, _, _ = np.linalg.lstsq(Ay, Y_arr[:, 2], rcond=None)
    else:
        m_h = 0.5 * (e_top[0] + e_bot[0])
        r_corr = Y_arr[:, 1] - m_h * Y_arr[:, 0]
        b1, a1 = np.polyfit(r_corr, Y_arr[:, 2], 1)
        ay = np.array([a1, b1, -b1 * m_h], float)
    uy_pred = ay[0] + ay[1] * Y_arr[:, 1] + ay[2] * Y_arr[:, 0]
    resid_y_px = float(np.max(np.abs((Y_arr[:, 2] - uy_pred) / ay[1])))

    return PanelCal(
        x_scale=x_scale,
        y_scale=y_scale,
        ax=ax,
        ay=ay,
        resid_x_px=resid_x_px,
        resid_y_px=resid_y_px,
        box=box,
        inner_box=inner,
    )


def mask_boxes(dark: np.ndarray, exclude_boxes_px: list[tuple[int, int, int, int]] | None) -> np.ndarray:
    """Return a copy of `dark` with excluded rectangular regions (text/legends) cleared."""
    if not exclude_boxes_px:
        return dark
    out = dark.copy()
    H, W = out.shape
    for c0, r0, c1, r1 in exclude_boxes_px:
        out[max(0, r0):min(H, r1), max(0, c0):min(W, c1)] = False
    return out


def detect_filled_markers(
    image: np.ndarray,
    cal: PanelCal,
    exclude_boxes_px: list[tuple[int, int, int, int]] | None = None,
    min_core_size_2x: int = 8,
    erode_2x: int = 7,
    min_edt_px: float = 6.5,
) -> list[tuple[float, float, int]]:
    """Fast filled-marker detector via 2x separable erosion + full-res EDT centroid refinement."""
    dark = mask_boxes(_to_dark(image), exclude_boxes_px)
    c0, r0, c1, r1 = cal.inner_box
    pad = 8
    roi = dark[r0 + pad:r1 - pad, c0 + pad:c1 - pad]
    roi2 = roi[::2, ::2]
    er = ndimage.binary_erosion(roi2, structure=np.ones((erode_2x, 1)))
    er = ndimage.binary_erosion(er, structure=np.ones((1, erode_2x)))
    lab, _ = ndimage.label(er)
    out = []
    for i, sl in enumerate(ndimage.find_objects(lab), start=1):
        if sl is None:
            continue
        sub = lab[sl] == i
        if int(sub.sum()) < min_core_size_2x:
            continue
        rr, cc = np.nonzero(sub)
        r_est = (rr.mean() + sl[0].start) * 2.0 + r0 + pad
        c_est = (cc.mean() + sl[1].start) * 2.0 + c0 + pad
        w = 34
        y0, y1 = max(r0 + 2, int(r_est) - w), min(r1 - 2, int(r_est) + w + 1)
        x0, x1 = max(c0 + 2, int(c_est) - w), min(c1 - 2, int(c_est) + w + 1)
        patch = dark[y0:y1, x0:x1]
        edt = ndimage.distance_transform_edt(patch)
        mx = float(edt.max())
        if mx < min_edt_px:
            continue
        core = edt >= max(min_edt_px - 1.0, 0.55 * mx)
        clab, cn = ndimage.label(core)
        best_d, best_pt = 1e9, None
        cy, cx = r_est - y0, c_est - x0
        for j in range(1, cn + 1):
            pr, pc = np.nonzero(clab == j)
            d2 = (pr.mean() - cy) ** 2 + (pc.mean() - cx) ** 2
            if d2 < best_d:
                best_d = d2
                best_pt = (float(pc.mean()) + x0, float(pr.mean()) + y0, int(len(pr)))
        if best_pt is not None:
            out.append(best_pt)
    return sorted(out, key=lambda p: p[0])


def detect_hollow_markers(
    image: np.ndarray,
    cal: PanelCal,
    exclude_boxes_px: list[tuple[int, int, int, int]] | None = None,
    fuse_px_2x: int = 5,
    min_hole_area_2x: int = 12,
    max_hole_area_2x: int = 800,
    min_cluster_area_2x: int = 20,
    max_bbox_2x: int = 44,
    pad: int = 6,
) -> list[tuple[float, float, int]]:
    """Fast hollow-marker detector via 2x enclosed white cavities fused across bisecting lines."""
    dark = mask_boxes(_to_dark(image), exclude_boxes_px)
    c0, r0, c1, r1 = cal.inner_box
    roi = dark[r0 + pad:r1 - pad, c0 + pad:c1 - pad]
    roi2 = roi[::2, ::2]

    wlab, _ = ndimage.label(~roi2)
    sizes = np.bincount(wlab.ravel())
    objs = ndimage.find_objects(wlab)
    small = np.zeros_like(roi2, dtype=bool)
    H2, W2 = roi2.shape
    for i, sl in enumerate(objs, start=1):
        if sl is None:
            continue
        a = int(sizes[i])
        if min_hole_area_2x <= a <= max_hole_area_2x:
            h = sl[0].stop - sl[0].start
            w = sl[1].stop - sl[1].start
            if h <= max_bbox_2x and w <= max_bbox_2x and sl[0].start > 0 and sl[0].stop < H2 and sl[1].start > 0 and sl[1].stop < W2:
                small[sl] |= (wlab[sl] == i)

    k = 2 * fuse_px_2x + 1
    fused = ndimage.binary_dilation(small, structure=np.ones((k, 1)))
    fused = ndimage.binary_dilation(fused, structure=np.ones((1, k)))
    hlab, _ = ndimage.label(fused)
    out = []
    for i, sl in enumerate(ndimage.find_objects(hlab), start=1):
        if sl is None:
            continue
        sub = small[sl] & (hlab[sl] == i)
        rr, cc = np.nonzero(sub)
        if len(rr) < min_cluster_area_2x:
            continue
        h = rr.max() - rr.min() + 1
        w = cc.max() - cc.min() + 1
        if not (6 <= h <= max_bbox_2x + 4 and 6 <= w <= max_bbox_2x + 4):
            continue
        r_est = (rr.mean() + sl[0].start) * 2.0 + r0 + pad
        c_est = (cc.mean() + sl[1].start) * 2.0 + c0 + pad
        hw = 38
        y0, y1 = max(r0 + 1, int(r_est) - hw), min(r1 - 1, int(r_est) + hw + 1)
        x0, x1 = max(c0 + 1, int(c_est) - hw), min(c1 - 1, int(c_est) + hw + 1)
        patch = dark[y0:y1, x0:x1]
        p_wlab, _ = ndimage.label(~patch)
        p_sizes = np.bincount(p_wlab.ravel())
        p_objs = ndimage.find_objects(p_wlab)
        pr_all, pc_all = [], []
        ph, pw = patch.shape
        cy, cx = r_est - y0, c_est - x0
        for j, psl in enumerate(p_objs, start=1):
            if psl is None:
                continue
            pa = int(p_sizes[j])
            if 20 <= pa <= 3200 and psl[0].start > 0 and psl[0].stop < ph and psl[1].start > 0 and psl[1].stop < pw:
                r_j, c_j = np.nonzero(p_wlab[psl] == j)
                r_j = r_j + psl[0].start
                c_j = c_j + psl[1].start
                if abs(r_j.mean() - cy) <= 22 and abs(c_j.mean() - cx) <= 22:
                    pr_all.append(r_j)
                    pc_all.append(c_j)
        if pr_all:
            rr_f = np.concatenate(pr_all)
            cc_f = np.concatenate(pc_all)
            c_ref = 0.5 * (0.5 * (cc_f.min() + cc_f.max()) + cc_f.mean()) + x0
            r_ref = 0.5 * (0.5 * (rr_f.min() + rr_f.max()) + rr_f.mean()) + y0
            out.append((float(c_ref), float(r_ref), int(len(rr_f))))
        else:
            out.append((float(c_est), float(r_est), int(len(rr) * 4)))
    return sorted(out, key=lambda p: p[0])


def segment_ink_fraction(dark: np.ndarray, p1: tuple[float, float], p2: tuple[float, float], margin_px: float = 26.0, perp_px: int = 3) -> float:
    """Fraction of sample points along straight line p1->p2 (outside margin_px of endpoints)
    that land within perp_px of dark ink."""
    c1, r1 = p1
    c2, r2 = p2
    dist = math.hypot(c2 - c1, r2 - r1)
    if dist <= 2.0 * margin_px + 4.0:
        return 1.0  # markers nearly touching
    n = max(12, int((dist - 2.0 * margin_px) / 4.0))
    ts = np.linspace(margin_px / dist, 1.0 - margin_px / dist, n)
    cs = np.rint(c1 + ts * (c2 - c1)).astype(int)
    rs = np.rint(r1 + ts * (r2 - r1)).astype(int)
    H, W = dark.shape
    hits = 0
    for c, r in zip(cs, rs):
        r0, r1_ = max(0, r - perp_px), min(H, r + perp_px + 1)
        c0, c1_ = max(0, c - perp_px), min(W, c + perp_px + 1)
        if dark[r0:r1_, c0:c1_].any():
            hits += 1
    return hits / n


def read_line_at(
    image: np.ndarray,
    cal: PanelCal,
    at: list[float],
    y_band: tuple[float, float] | None = None,
    exclude_boxes_px: list[tuple[int, int, int, int]] | None = None,
    half_win_c: int = 3,
    min_thick: int = 2,
    max_thick: int = 18,
    note_prefix: str = "fit line",
) -> list[dict]:
    """Read a thin straight or smooth black line at explicit x values `at`."""
    dark = mask_boxes(_to_dark(image), exclude_boxes_px)
    c0, r0, c1, r1 = cal.inner_box
    out = []
    for xv in at:
        c_mid = cal.col_for_x(xv, 0.5 * (r0 + r1))
        if y_band is not None:
            ra = cal.row_for_y(y_band[0], c_mid)
            rb = cal.row_for_y(y_band[1], c_mid)
            r_lo = max(r0 + 12, int(round(min(ra, rb))))
            r_hi = min(r1 - 12, int(round(max(ra, rb))))
        else:
            r_lo, r_hi = r0 + 12, r1 - 12

        rows_found = []
        thicks = []
        for dc in range(-half_win_c, half_win_c + 1):
            cc = int(round(c_mid)) + dc
            if cc <= c0 + 4 or cc >= c1 - 4:
                continue
            col_dark = dark[r_lo:r_hi, cc]
            idx = np.nonzero(col_dark)[0]
            if len(idx) == 0:
                continue
            for gr in np.split(idx, np.where(np.diff(idx) > 2)[0] + 1):
                th = int(gr[-1] - gr[0] + 1)
                if min_thick <= th <= max_thick:
                    rows_found.append(r_lo + 0.5 * (gr[0] + gr[-1]))
                    thicks.append(th)
        if not rows_found:
            out.append({"x": float(xv), "y": None, "note": f"{note_prefix}: line not found"})
            continue
        r_med = float(np.median(rows_found))
        c_exact = cal.col_for_x(xv, r_med)
        half_px = max(2.5, 0.5 * float(np.median(thicks)) + 1.0 + cal.resid_y_px)
        y, lo, hi, px_lo, px_hi = cal.y_interval(c_exact, r_med, half_px)
        out.append({
            "x": float(xv),
            "y": y,
            "lo": lo,
            "hi": hi,
            "px": [c_exact, r_med],
            "px_lo": px_lo,
            "px_hi": px_hi,
            "note": f"{note_prefix} (thick {np.median(thicks):.1f} px, cal resid {cal.resid_y_px:.2f} px)",
        })
    return out


def make_geom_template(shape: str, radius_px: float = 28.0, stroke_px: float = 7.0) -> tuple[np.ndarray, np.ndarray]:
    """Build exact full-resolution geometric template (pos_ink, neg_white) for KaleidaGraph symbols."""
    R = float(radius_px)
    S = int(round(2 * R)) + 15
    if S % 2 == 0:
        S += 1
    cy = cx = S // 2
    yy, xx = np.mgrid[0:S, 0:S] - cy
    rr = np.hypot(xx, yy)
    t = stroke_px * 0.6

    if shape == "circle":
        pos = np.abs(rr - R) <= t
        neg = rr < (R - stroke_px * 0.9)
        return pos, neg
    if shape == "dotcircle":
        pos = (np.abs(rr - R) <= t) | (rr <= 5.5)
        neg = (rr < (R - stroke_px * 0.9)) & (rr > 8.0)
        return pos, neg
    if shape == "square":
        m = np.maximum(np.abs(xx), np.abs(yy))
        pos = np.abs(m - R * 0.92) <= t
        neg = m < (R * 0.92 - stroke_px * 0.9)
        return pos, neg
    if shape == "diamond":
        m = np.abs(xx) + np.abs(yy)
        pos = np.abs(m - R * 1.15) <= t * 1.25
        neg = m < (R * 1.15 - stroke_px * 1.2)
        return pos, neg
    if shape in ("uptri", "downtri"):
        # KaleidaGraph triangle: flat base width ~2*R, total height ~2*R, center at midpoint of height
        y_top, y_bot = -R, R
        half_w = (yy - y_top) / (y_bot - y_top) * R
        inside = (yy >= y_top) & (yy <= y_bot) & (np.abs(xx) <= half_w)
        d_bot = np.abs(yy - y_bot)
        slope = R / (y_bot - y_top)
        d_side = np.abs(np.abs(xx) - half_w) / math.hypot(1.0, slope)
        pos = ((d_bot <= t) & (np.abs(xx) <= R + t)) | ((d_side <= t) & (yy >= y_top - t) & (yy <= y_bot + t))
        neg = inside & (d_bot > stroke_px * 0.9) & (d_side > stroke_px * 0.9)
        if shape == "downtri":
            return pos[::-1, :].copy(), neg[::-1, :].copy()
        return pos, neg
    if shape == "filled_circle":
        pos = rr <= R * 0.9
        neg = (rr >= R * 1.2) & (rr <= R * 1.55)
        return pos, neg
    if shape in ("filled_uptri", "filled_downtri"):
        y_top, y_bot = -R, R
        half_w = (yy - y_top) / (y_bot - y_top) * R
        pos = (yy >= y_top) & (yy <= y_bot) & (np.abs(xx) <= half_w)
        neg = (~pos) & (yy >= y_top - 8) & (yy <= y_bot + 8) & (np.abs(xx) <= R + 8)
        if shape == "filled_downtri":
            return pos[::-1, :].copy(), neg[::-1, :].copy()
        return pos, neg
    raise ValueError(f"Unknown template shape: {shape}")


def score_patch_template(dark: np.ndarray, c: float, r: float, pos: np.ndarray, neg: np.ndarray, neg_w: float = 0.4, search_r: int = 8) -> tuple[float, float, float]:
    """Local correlation search around (c, r) within ±search_r px (step 2 px) using exact full-res template.
    Returns (best_score, best_c, best_r)."""
    H, W = dark.shape
    sh, sw = pos.shape
    hy, hx = sh // 2, sw // 2
    ci, ri = int(round(c)), int(round(r))
    y0, y1 = max(0, ri - hy - search_r), min(H, ri + hy + search_r + 1)
    x0, x1 = max(0, ci - hx - search_r), min(W, ci + hx + search_r + 1)
    sub = dark[y0:y1, x0:x1].astype(np.float32)
    if sub.shape[0] < sh or sub.shape[1] < sw:
        return -1.0, c, r
    k = pos.astype(np.float32) / max(1, int(pos.sum())) - neg_w * neg.astype(np.float32) / max(1, int(neg.sum()))
    resp = signal.fftconvolve(sub, k[::-1, ::-1], mode="valid")
    idx = np.unravel_index(int(np.argmax(resp)), resp.shape)
    best_r = float(y0 + idx[0] + hy)
    best_c = float(x0 + idx[1] + hx)
    return float(resp[idx]), best_c, best_r


def find_template_in_band(
    dark: np.ndarray,
    c_range: tuple[int, int],
    r_range: tuple[int, int],
    pos: np.ndarray,
    neg: np.ndarray,
    neg_w: float = 0.25,
    step: int = 2,
) -> tuple[float, float, float]:
    """Find best match of full-res template inside a pixel bounding box (c0..c1, r0..r1).
    Evaluates at `step`x downsample for speed, then refines locally at 1 px."""
    H, W = dark.shape
    c0, c1 = max(0, c_range[0]), min(W, c_range[1])
    r0, r1 = max(0, r_range[0]), min(H, r_range[1])
    sub = dark[r0:r1:step, c0:c1:step].astype(np.float32)
    p2 = pos[::step, ::step]
    n2 = neg[::step, ::step]
    if sub.shape[0] < p2.shape[0] or sub.shape[1] < p2.shape[1]:
        return -1.0, 0.5 * (c0 + c1), 0.5 * (r0 + r1)
    k = p2.astype(np.float32) / max(1, int(p2.sum())) - neg_w * n2.astype(np.float32) / max(1, int(n2.sum()))
    resp = signal.fftconvolve(sub, k[::-1, ::-1], mode="valid")
    idx = np.unravel_index(int(np.argmax(resp)), resp.shape)
    r_est = r0 + (idx[0] + p2.shape[0] // 2) * step
    c_est = c0 + (idx[1] + p2.shape[1] // 2) * step
    return score_patch_template(dark, c_est, r_est, pos, neg, neg_w=neg_w, search_r=4)


def dilate_for_trace(dark: np.ndarray, radius: int = 4) -> np.ndarray:
    """Fast separable dilation of dark mask for vectorized line-segment ink checks."""
    k = 2 * radius + 1
    out = ndimage.binary_dilation(dark, structure=np.ones((k, 1), dtype=bool))
    return ndimage.binary_dilation(out, structure=np.ones((1, k), dtype=bool))


def fast_segment_ink(dark_dil: np.ndarray, p1: tuple[float, float], p2: tuple[float, float], margin_px: float = 34.0) -> float:
    """Vectorized fraction of samples along straight segment p1->p2 hitting dark_dil."""
    c1, r1 = p1
    c2, r2 = p2
    dist = math.hypot(c2 - c1, r2 - r1)
    if dist <= 2.0 * margin_px + 4.0:
        return 1.0
    n = max(12, int((dist - 2.0 * margin_px) / 4.0))
    ts = np.linspace(margin_px / dist, 1.0 - margin_px / dist, n)
    cs = np.clip(np.rint(c1 + ts * (c2 - c1)).astype(int), 0, dark_dil.shape[1] - 1)
    rs = np.clip(np.rint(r1 + ts * (r2 - r1)).astype(int), 0, dark_dil.shape[0] - 1)
    return float(dark_dil[rs, cs].mean())


def column_clusters(dark: np.ndarray, cal: PanelCal, x: float, y_band: tuple[float, float] | None = None,
                    half_c: int = 3, max_gap: int = 4) -> list[tuple[float, int, int]]:
    """Return list of (r_center, r_top, r_bot) dark runs along the slanted column for x."""
    c0, r0, c1, r1 = cal.inner_box
    c_mid = cal.col_for_x(x, 0.5 * (r0 + r1))
    if y_band is not None:
        ra, rb = cal.row_for_y(y_band[0], c_mid), cal.row_for_y(y_band[1], c_mid)
        r_lo, r_hi = max(r0 + 8, int(round(min(ra, rb)))), min(r1 - 8, int(round(max(ra, rb))))
    else:
        r_lo, r_hi = r0 + 8, r1 - 8
    ci = int(round(c_mid))
    strip = dark[r_lo:r_hi, max(c0 + 2, ci - half_c):min(c1 - 2, ci + half_c + 1)]
    hits = np.nonzero(strip.any(axis=1))[0]
    if len(hits) == 0:
        return []
    out = []
    for g in np.split(hits, np.where(np.diff(hits) > max_gap)[0] + 1):
        out.append((r_lo + 0.5 * (g[0] + g[-1]), r_lo + int(g[0]), r_lo + int(g[-1])))
    return out


def row_clusters(dark: np.ndarray, cal: PanelCal, y: float, x_band: tuple[float, float] | None = None,
                 half_r: int = 3, max_gap: int = 4) -> list[tuple[float, int, int]]:
    """Return list of (c_center, c_left, c_right) dark runs along the horizontal row for y."""
    c0, r0, c1, r1 = cal.inner_box
    r_mid = cal.row_for_y(y, 0.5 * (c0 + c1))
    if x_band is not None:
        ca, cb = cal.col_for_x(x_band[0], r_mid), cal.col_for_x(x_band[1], r_mid)
        c_lo, c_hi = max(c0 + 8, int(round(min(ca, cb)))), min(c1 - 8, int(round(max(ca, cb))))
    else:
        c_lo, c_hi = c0 + 8, c1 - 8
    ri = int(round(r_mid))
    strip = dark[max(r0 + 2, ri - half_r):min(r1 - 2, ri + half_r + 1), c_lo:c_hi]
    hits = np.nonzero(strip.any(axis=0))[0]
    if len(hits) == 0:
        return []
    out = []
    for g in np.split(hits, np.where(np.diff(hits) > max_gap)[0] + 1):
        out.append((c_lo + 0.5 * (g[0] + g[-1]), c_lo + int(g[0]), c_lo + int(g[-1])))
    return out


def trace_step_along_curve(
    dark: np.ndarray,
    cal: PanelCal,
    c_start: float,
    r_start: float,
    c_target: float,
    r_pred_slope: float = 0.0,
    search_r_px: int = 55,
    step_c: float = 18.0,
) -> tuple[float, float]:
    """Walk column-by-column along a continuous dark curve from (c_start, r_start) to c_target,
    maintaining a smoothed slope so crossings and hollow markers are traversed cleanly.
    Returns (c_target, r_end)."""
    c0, r0, c1, r1 = cal.inner_box
    direction = 1.0 if c_target >= c_start else -1.0
    c = float(c_start)
    r = float(r_start)
    slope = float(r_pred_slope)  # dr/dc
    n_steps = max(1, int(math.ceil(abs(c_target - c_start) / step_c)))
    dc = (c_target - c_start) / n_steps
    hist_c, hist_r = [c], [r]
    for _ in range(n_steps):
        c_next = c + dc
        r_pred = r + slope * dc
        ci = int(round(c_next))
        if ci <= c0 + 2 or ci >= c1 - 2:
            c, r = c_next, r_pred
            continue
        r_lo = max(r0 + 4, int(round(r_pred)) - search_r_px)
        r_hi = min(r1 - 4, int(round(r_pred)) + search_r_px + 1)
        strip = dark[r_lo:r_hi, max(c0 + 2, ci - 2):min(c1 - 2, ci + 3)]
        hits = np.nonzero(strip.any(axis=1))[0]
        if len(hits) == 0:
            c, r = c_next, r_pred
            continue
        # Group runs; also merge symmetric top/bottom rims of a hollow marker if straddling r_pred
        groups = [r_lo + 0.5 * (g[0] + g[-1]) for g in np.split(hits, np.where(np.diff(hits) > 4)[0] + 1)]
        # If two runs bracket r_pred within 20..65 px (hollow marker!), take their midpoint!
        r_best = None
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                sep = groups[j] - groups[i]
                mid = 0.5 * (groups[i] + groups[j])
                if 20.0 <= sep <= 72.0 and abs(mid - r_pred) <= 22.0:
                    if r_best is None or abs(mid - r_pred) < abs(r_best - r_pred):
                        r_best = mid
        if r_best is None:
            r_best = min(groups, key=lambda rm: abs(rm - r_pred))
        # Update smoothed slope over recent points
        hist_c.append(c_next)
        hist_r.append(r_best)
        if len(hist_c) >= 3:
            k_win = min(5, len(hist_c))
            dc_w = hist_c[-1] - hist_c[-k_win]
            if abs(dc_w) > 10.0:
                slope = (hist_r[-1] - hist_r[-k_win]) / dc_w
                slope = max(-6.0, min(6.0, slope))
        c, r = c_next, r_best
    return float(c), float(r)


def refine_hollow_or_ink_at(
    dark: np.ndarray,
    cal: PanelCal,
    c_est: float,
    r_est: float,
    win: int = 38,
) -> tuple[float, float, float]:
    """Refine a marker center near (c_est, r_est) from enclosed white cavities (if clean hollow
    marker) or from vertical ink symmetry around (c_est, r_est). Returns (c, r, half_err_px)."""
    c0, r0, c1, r1 = cal.inner_box
    y0, y1 = max(r0 + 2, int(round(r_est)) - win), min(r1 - 2, int(round(r_est)) + win + 1)
    x0, x1 = max(c0 + 2, int(round(c_est)) - win), min(c1 - 2, int(round(c_est)) + win + 1)
    patch = dark[y0:y1, x0:x1]
    ph, pw = patch.shape
    cy, cx = r_est - y0, c_est - x0
    # Check for enclosed white cavities within 20 px of (cx, cy)
    p_d = ndimage.binary_dilation(patch, iterations=1)
    wlab, _ = ndimage.label(~p_d)
    sizes = np.bincount(wlab.ravel())
    objs = ndimage.find_objects(wlab)
    pr_all, pc_all = [], []
    for j, sl in enumerate(objs, start=1):
        if sl is None:
            continue
        a = int(sizes[j])
        if 20 <= a <= 3000 and sl[0].start > 0 and sl[0].stop < ph and sl[1].start > 0 and sl[1].stop < pw:
            h = sl[0].stop - sl[0].start
            w = sl[1].stop - sl[1].start
            if h <= 65 and w <= 65:
                rr, cc = np.nonzero(wlab[sl] == j)
                rr = rr + sl[0].start
                cc = cc + sl[1].start
                if math.hypot(rr.mean() - cy, cc.mean() - cx) <= 22.0:
                    pr_all.append(rr)
                    pc_all.append(cc)
    if pr_all:
        rr_f = np.concatenate(pr_all)
        cc_f = np.concatenate(pc_all)
        c_ref = 0.5 * (0.5 * (cc_f.min() + cc_f.max()) + cc_f.mean()) + x0
        r_ref = 0.5 * (0.5 * (rr_f.min() + rr_f.max()) + rr_f.mean()) + y0
        if math.hypot(c_ref - c_est, r_ref - r_est) <= 20.0:
            return float(c_ref), float(r_ref), max(3.5, 2.0 + cal.resid_y_px)
    return float(c_est), float(r_est), max(4.5, 3.0 + cal.resid_y_px)

