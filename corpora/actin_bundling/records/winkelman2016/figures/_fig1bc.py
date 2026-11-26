"""Calibrated reader for Fig. 1B (fascin) and Fig. 1C (alpha-actinin) low-speed sedimentation."""

from __future__ import annotations

import numpy as np
from scipy import ndimage

import _calib as C


def _calibrate_1b(image: np.ndarray):
    dark, gray, h_lines, v_lines = C.find_l_spines(image, thresh=140, row_frac=0.30, col_frac=0.30)
    x_spine = max(h_lines, key=lambda t: t[0])
    y_spine = min(v_lines, key=lambda t: t[0])

    yt_raw = C.find_ticks(dark, y_spine, (100, int(x_spine[1]) - 40), "y", side=+1, band=(4, 24), frac=0.4, max_width=40)
    yt = C.match_ticks_regular(yt_raw, 3, tol=6.0, min_step=300.0)
    y_axis = C.verify_axis("linear", list(zip(yt, [0.15, 0.10, 0.05])), max_residual_px=2.0)

    xt_raw = C.find_ticks(dark, x_spine, (int(y_spine[2]) + 100, image.shape[1] - 100), "x", side=-1, band=(4, 24), frac=0.4, max_width=40)
    xt = C.match_ticks_regular(xt_raw, 5, tol=8.0, min_step=250.0)
    x_axis = C.verify_axis("linear", list(zip(xt, [0.2, 0.4, 0.6, 0.8, 1.0])), max_residual_px=3.0)
    return x_axis, y_axis, x_spine, y_spine


def _calibrate_1c(image: np.ndarray):
    dark, gray, h_lines, v_lines = C.find_l_spines(image, thresh=140, row_frac=0.30, col_frac=0.30)
    x_spine = max(h_lines, key=lambda t: t[0])
    y_spine = min(v_lines, key=lambda t: t[0])

    yt_raw = C.find_ticks(dark, y_spine, (40, int(x_spine[1]) - 30), "y", side=+1, band=(4, 24), frac=0.4, max_width=40)
    yt = C.match_ticks_regular(yt_raw, 5, tol=6.0, min_step=300.0)
    y_axis = C.verify_axis("linear", list(zip(yt, [0.08, 0.06, 0.04, 0.02, 0.0])), max_residual_px=2.0)

    xt_raw = C.find_ticks(dark, x_spine, (int(y_spine[2]) + 20, image.shape[1] - 20), "x", side=-1, band=(4, 24), frac=0.4, max_width=40)
    xt = C.match_ticks_regular(xt_raw, 7, tol=8.0, min_step=250.0)
    x_axis = C.verify_axis("linear", list(zip(xt, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])), max_residual_px=3.0)
    return x_axis, y_axis, x_spine, y_spine


def _find_circle_centers_fast(mask: np.ndarray, min_dist_px: float = 24.0, min_peak_px: float = 21.0) -> list[tuple[float, float]]:
    """Fast circle center finder: 2x-downsampled distance transform + full-res local disk refinement."""
    sub2 = mask[::2, ::2]
    dist2 = ndimage.distance_transform_edt(sub2)
    sz = max(3, int(round(min_dist_px)) | 1)
    loc_max = (dist2 == ndimage.maximum_filter(dist2, size=sz)) & (dist2 >= min_peak_px / 2.0)
    lbl, n = ndimage.label(loc_max)
    pts = []
    H, W = mask.shape
    for i in range(1, n + 1):
        rr, cc = np.where(lbl == i)
        r0, c0 = float(rr.mean()) * 2.0, float(cc.mean()) * 2.0
        rad = int(round(dist2[int(round(rr.mean())), int(round(cc.mean()))] * 2.0))
        y0, y1 = max(0, int(r0) - rad), min(H, int(r0) + rad + 1)
        x0, x1 = max(0, int(c0) - rad), min(W, int(c0) + rad + 1)
        yy, xx = np.where(mask[y0:y1, x0:x1])
        d2 = (yy + y0 - r0) ** 2 + (xx + x0 - c0) ** 2
        keep = d2 <= (rad * 0.90) ** 2
        if keep.sum() > 20:
            pts.append((float((xx[keep] + x0).mean()) + 0.5, float((yy[keep] + y0).mean()) + 0.5))
        else:
            pts.append((c0 + 0.5, r0 + 0.5))
    pts.sort(key=lambda p: p[0])
    return pts


def _read_markers(image: np.ndarray, panel: str) -> list[dict]:
    x_axis, y_axis, x_spine, y_spine = _calibrate_1b(image) if panel == "1b" else _calibrate_1c(image)
    R = image[:, :, 0].astype(np.float32)
    G = image[:, :, 1].astype(np.float32)
    B = image[:, :, 2].astype(np.float32)
    if panel == "1b":
        cmask = (R > 120) & (R - G > 50) & (R - B > 50)
        row_min = 60
        rad_mark = 33.0
    else:
        cmask = (B > 120) & (B - R > 50) & (B - G > 20)
        row_min = 40
        rad_mark = 36.0
    row_max = int(x_spine[1]) - 2
    col_min = int(y_spine[2]) + 5
    cmask[:row_min, :] = False
    cmask[row_max + 1 :, :] = False
    cmask[:, :col_min] = False

    centers = _find_circle_centers_fast(cmask, min_dist_px=24.0, min_peak_px=21.0)
    # Build a mask of the markers themselves so error-bar tracing does not walk into an adjacent dot
    marker_only = np.zeros_like(cmask)
    yy_grid, xx_grid = np.ogrid[: cmask.shape[0], : cmask.shape[1]]
    for c_m, r_m in centers:
        y0, y1 = max(0, int(r_m) - 42), min(cmask.shape[0], int(r_m) + 43)
        x0, x1 = max(0, int(c_m) - 42), min(cmask.shape[1], int(c_m) + 43)
        sub_d2 = (yy_grid[y0:y1, :] - r_m) ** 2 + (xx_grid[:, x0:x1] - c_m) ** 2
        marker_only[y0:y1, x0:x1] |= sub_d2 <= (rad_mark + 4.0) ** 2

    out = []
    for idx, (c_mid, r_mid) in enumerate(centers):
        x_val = float(x_axis.to_value(c_mid))
        if abs(x_val) < 0.005:
            x_val = 0.0
        # Exclude OTHER markers from stem_mask so closely spaced dots do not falsely extend bars
        stem_mask = cmask.copy()
        for j, (cj, rj) in enumerate(centers):
            if j == idx:
                continue
            y0, y1 = max(0, int(rj) - 42), min(cmask.shape[0], int(rj) + 43)
            x0, x1 = max(0, int(cj) - 42), min(cmask.shape[1], int(cj) + 43)
            sub_d2 = (yy_grid[y0:y1, :] - rj) ** 2 + (xx_grid[:, x0:x1] - cj) ** 2
            stem_mask[y0:y1, x0:x1] &= ~(sub_d2 <= (rad_mark + 2.0) ** 2)

        u_row, l_row = C.find_vertical_error_bar(
            stem_mask,
            c_mid,
            r_mid,
            row_min=row_min,
            row_max=row_max,
            stem_half_w=8,
            max_gap=8,
            cap_half_w=45,
            min_cap_width=30,
            start_up_offset=int(rad_mark - 4),
            start_dn_offset=int(rad_mark - 4),
        )
        if u_row is not None and (r_mid - u_row) < rad_mark + 8:
            u_row = None
        if l_row is not None and (l_row - r_mid) < rad_mark + 8:
            l_row = None
        # Symmetrize if only one side found (when other side was occluded by neighboring dot or x-spine)
        if u_row is not None and l_row is None:
            l_row = r_mid + (r_mid - u_row)
        elif l_row is not None and u_row is None:
            u_row = r_mid - (l_row - r_mid)

        note = f"Fig. {panel.upper()} marker; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px"
        if u_row is None and l_row is None:
            rd = C.make_point_reading(
                x_val,
                c_mid,
                r_mid,
                y_axis,
                upper_err_row=r_mid - rad_mark,
                lower_err_row=r_mid + rad_mark,
                mirror_lower=False,
                err_px=2.0 + y_axis.residual_px,
                note=note + "; SEM bar within marker radius",
            )
        else:
            rd = C.make_point_reading(
                x_val,
                c_mid,
                r_mid,
                y_axis,
                upper_err_row=u_row,
                lower_err_row=l_row,
                mirror_lower=False,
                err_px=2.0 + y_axis.residual_px,
                note=note,
            )
        out.append(rd)
    return out


def _read_fit(image: np.ndarray, panel: str, at: list[float]) -> list[dict]:
    x_axis, y_axis, x_spine, y_spine = _calibrate_1b(image) if panel == "1b" else _calibrate_1c(image)
    R = image[:, :, 0].astype(np.float32)
    G = image[:, :, 1].astype(np.float32)
    B = image[:, :, 2].astype(np.float32)
    blk = (R < 90) & (G < 90) & (B < 90)
    row_min = 60
    row_max = int(x_spine[1]) - 15
    col_min = int(y_spine[2]) + 15
    blk[:row_min, :] = False
    blk[row_max + 1 :, :] = False
    blk[:, :col_min] = False

    # Trace curve right-to-left skipping marker occlusions cleanly
    xs_dense = np.linspace(1.18, 0.08, 160)
    trace_c, trace_r = [], []
    prev_r = None
    for xv in xs_dense:
        c = int(round(x_axis.to_pixel(xv)))
        if c < col_min or c >= image.shape[1]:
            continue
        rs = C.runs_1d(blk[row_min : row_max + 1, c])
        rs = [(a + row_min, b + row_min) for a, b in rs if (b - a + 1) <= 55]
        if not rs:
            continue
        if prev_r is None:
            seg = min(rs, key=lambda t: t[0])
        else:
            # Expect monotonically increasing row (descending curve from right to left)
            cands = [s for s in rs if (0.5 * (s[0] + s[1]) - prev_r) >= -20 and (0.5 * (s[0] + s[1]) - prev_r) <= 220]
            if not cands:
                continue
            seg = min(cands, key=lambda t: abs(0.5 * (t[0] + t[1]) - prev_r))
        prev_r = 0.5 * (seg[0] + seg[1]) + 0.5
        trace_c.append(float(c) + 0.5)
        trace_r.append(prev_r)
    trace_c = np.array(trace_c[::-1])
    trace_r = np.array(trace_r[::-1])

    out = []
    for xv in at:
        c_f = float(x_axis.to_pixel(xv))
        pred_r = float(np.interp(c_f, trace_c, trace_r))

        def pick(rs, col, _pr=pred_r):
            cands = [s for s in rs if (s[1] - s[0] + 1) <= 55 and abs(0.5 * (s[0] + s[1]) + 0.5 - _pr) <= 75]
            if not cands:
                return None
            return min(cands, key=lambda s: abs(0.5 * (s[0] + s[1]) + 0.5 - _pr))

        rd = C.sample_solid_curve_y(
            blk,
            x_axis,
            y_axis,
            float(xv),
            row_min,
            row_max,
            pick_fn=pick,
            search_dx=10,
            note_prefix=f"Fig. {panel.upper()} Hill fit; ",
        )
        out.append(rd)
    return out


def read_series(image: np.ndarray, panel: str, kind: str, at: list[float]) -> list[dict]:
    if kind == "data":
        return _read_markers(image, panel)
    if kind == "fit":
        return _read_fit(image, panel, at)
    raise ValueError(f"Unknown kind {kind!r}")
