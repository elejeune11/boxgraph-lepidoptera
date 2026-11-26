"""Calibrated reader for Fig. 3H (fascin domain length vs [Fascin]/[alpha-Actinin]: data, coop & simple models)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage

import _calib as C


def _calibrate(image: np.ndarray):
    dark, gray, h_lines, v_lines = C.find_l_spines(
        image, thresh=100, row_frac=0.35, col_frac=0.35, row_range=(120, 1800), col_range=(350, 2200)
    )
    x_spine = max(h_lines, key=lambda t: t[0])
    y_spine = min(v_lines, key=lambda t: t[0])

    # Inward y ticks: 8, 6, 4, 2 (0 is overlapped by red data point at (0,0))
    yt_raw = C.find_ticks(dark, y_spine, (180, int(x_spine[1]) - 150), "y", side=+1, band=(4, 24), frac=0.35, max_width=45)
    yt = C.match_ticks_regular(yt_raw, 4, tol=8.0, min_step=250.0)
    y_axis = C.verify_axis("linear", list(zip(yt, [8.0, 6.0, 4.0, 2.0])), max_residual_px=2.5)

    # Inward x ticks: 0, 0.4, 0.8, 1.2, 1.6
    xt_raw = C.find_ticks(dark, x_spine, (int(y_spine[2]) + 15, image.shape[1] - 30), "x", side=-1, band=(4, 24), frac=0.35, max_width=45)
    xt = C.match_ticks_regular(xt_raw, 5, tol=8.0, min_step=300.0)
    x_axis = C.verify_axis("linear", list(zip(xt, [0.0, 0.4, 0.8, 1.2, 1.6])), max_residual_px=2.5)
    return x_axis, y_axis, x_spine, y_spine


def _find_legend_box(image: np.ndarray, y_spine) -> tuple[int, int, int, int]:
    """Find the black legend rectangle in upper-left of plot box."""
    dark = image.astype(np.float32).mean(axis=2) < 100
    c_lo = int(y_spine[2]) + 30
    c_hi = min(image.shape[1], c_lo + 1150)
    sub = dark[180:660, c_lo:c_hi]
    h_rows = np.where(sub.sum(axis=1) > 700)[0] + 180
    if len(h_rows) < 2:
        return (0, 0, 0, 0)
    r0, r1 = int(h_rows.min()), int(h_rows.max())
    v_cols = np.where(dark[r0 : r1 + 1, c_lo:c_hi].sum(axis=0) > 0.7 * (r1 - r0))[0] + c_lo
    if len(v_cols) < 2:
        return (0, 0, 0, 0)
    return (int(v_cols.min()) - 6, r0 - 6, int(v_cols.max()) + 6, r1 + 6)


def _read_data(image: np.ndarray) -> list[dict]:
    x_axis, y_axis, x_spine, y_spine = _calibrate(image)
    R = image[:, :, 0].astype(np.float32)
    G = image[:, :, 1].astype(np.float32)
    B = image[:, :, 2].astype(np.float32)
    col_min = int(y_spine[2]) + 5
    # Find where top micrograph black border ends (rows 0..140)
    top_dark = (image[:150, col_min:].astype(np.float32).mean(axis=2) < 200).mean(axis=1) > 0.2
    top_rows = np.where(top_dark)[0]
    row_min = int(top_rows.max()) + 5 if len(top_rows) > 0 else 95
    row_max = int(x_spine[1]) - 5

    red = (R > 160) & (R - G > 80) & (R - B > 80)
    red[:row_min, :] = False
    red[row_max + 1 :, :] = False
    red[:, :col_min] = False

    lc0, lr0, lc1, lr1 = _find_legend_box(image, y_spine)
    if lc1 > lc0 and lr1 > lr0:
        red[lr0:lr1, lc0:lc1] = False

    lbl, n = ndimage.label(red)
    pts = []
    for i in range(1, n + 1):
        rr, cc = np.where(lbl == i)
        if len(rr) > 500:
            pts.append((float(cc.mean()) + 0.5, float(rr.mean()) + 0.5))
    pts.sort(key=lambda p: p[0])

    # Black mask for vertical error bars (exclude legend box, spines, and top micrograph border)
    blk = (R < 85) & (G < 85) & (B < 85)
    eb_row_min = row_min
    blk[:eb_row_min, :] = False
    blk[int(x_spine[1]) - 2 :, :] = False
    blk[:, :col_min] = False
    if lc1 > lc0 and lr1 > lr0:
        blk[lr0:lr1, lc0:lc1] = False

    # Build pure-vertical stem mask: keep black pixels that lie in a vertical column of height >= 10
    # or horizontal cap of width >= 22 (to prevent walking along sloped connecting polylines!)
    rad_mark = 42.0
    out = []
    for c_mid, r_mid in pts:
        xv = float(x_axis.to_value(c_mid))
        if abs(xv) < 0.025:
            xv = 0.0
        # Trace strictly within narrow vertical strip around c_mid (stem width ~10 px)
        c0 = max(col_min, int(round(c_mid)) - 6)
        c1 = min(blk.shape[1], int(round(c_mid)) + 7)
        # Require at least 4 dark pixels across the 13-px strip (avoids thin angled polylines leaving c_mid)
        strip = blk[:, c0:c1].sum(axis=1) >= 4
        # Also allow bridging up to 10 px where blue curve or red dot crosses
        u_end = None
        gap = 0
        for r in range(max(eb_row_min, int(round(r_mid)) - 34), eb_row_min - 1, -1):
            if strip[r]:
                u_end = r
                gap = 0
            else:
                gap += 1
                if gap > 12:
                    break
        l_end = None
        gap = 0
        for r in range(min(row_max, int(round(r_mid)) + 34), row_max + 1):
            if strip[r]:
                l_end = r
                gap = 0
            else:
                gap += 1
                if gap > 12:
                    break

        def refine_cap(end_r, direction):
            if end_r is None:
                return None
            k0 = max(col_min, int(round(c_mid)) - 32)
            k1 = min(blk.shape[1], int(round(c_mid)) + 33)
            lo_r = max(eb_row_min, end_r - 10)
            hi_r = min(row_max, end_r + 10)
            widths = blk[lo_r : hi_r + 1, k0:k1].sum(axis=1)
            cap_rows = np.where(widths >= 24)[0] + lo_r
            if len(cap_rows) > 0:
                cs = C.clusters(cap_rows, gap=2)
                best = min(cs, key=lambda t: t[0]) if direction < 0 else max(cs, key=lambda t: t[0])
                return float(best[0])
            return float(end_r) + 0.5

        u_row = refine_cap(u_end, -1)
        l_row = refine_cap(l_end, +1)
        # Filter tiny spurious hits within marker radius (46 px)
        if u_row is not None and (r_mid - u_row) < 46:
            u_row = None
        if l_row is not None and (l_row - r_mid) < 46:
            l_row = None
        # Note: if point is at x=0 (y=0), no error bar is plotted
        if xv == 0.0:
            u_row, l_row = None, None
        elif u_row is not None and l_row is None:
            l_row = r_mid + (r_mid - u_row)
        elif l_row is not None and u_row is None:
            u_row = r_mid - (l_row - r_mid)

        note = f"Fig. 3H data; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px"
        if xv == 0.0:
            rd = C.make_point_reading(
                xv, c_mid, r_mid, y_axis, err_px=2.5 + y_axis.residual_px, note=note + "; zero-fascin control (no error bar)"
            )
        elif u_row is None and l_row is None:
            rd = C.make_point_reading(
                xv,
                c_mid,
                r_mid,
                y_axis,
                upper_err_row=r_mid - rad_mark,
                lower_err_row=r_mid + rad_mark,
                mirror_lower=False,
                err_px=2.5 + y_axis.residual_px,
                note=note + "; SEM bar within marker radius",
            )
        else:
            rd = C.make_point_reading(
                xv,
                c_mid,
                r_mid,
                y_axis,
                upper_err_row=u_row,
                lower_err_row=l_row,
                mirror_lower=False,
                err_px=2.5 + y_axis.residual_px,
                note=note,
            )
        out.append(rd)
    return out


def _read_model_curve(image: np.ndarray, mode: str, at: list[float]) -> list[dict]:
    x_axis, y_axis, x_spine, y_spine = _calibrate(image)
    R = image[:, :, 0].astype(np.float32)
    G = image[:, :, 1].astype(np.float32)
    B = image[:, :, 2].astype(np.float32)
    blue = (B > 150) & (B - R > 80) & (B - G > 50)
    row_min = 120
    row_max = int(x_spine[1]) - 5
    col_min = int(y_spine[2]) + 5
    blue[:row_min, :] = False
    blue[row_max + 1 :, :] = False
    blue[:, :col_min] = False

    lc0, lr0, lc1, lr1 = _find_legend_box(image, y_spine)
    if lc1 > lc0 and lr1 > lr0:
        blue[lr0:lr1, lc0:lc1] = False

    # Row threshold separating simple model (y < 0.35 um) from coop model (y > 0.35 um for x >= 0.1)
    r_split = int(round(y_axis.to_pixel(0.35)))

    if mode == "coop":
        # Solid upper blue curve
        coop_mask = blue.copy()
        coop_mask[r_split + 40 :, :] = False
        out = []
        for xv in at:
            r_hi = row_max if xv < 0.15 else r_split + 20
            rd = C.sample_solid_curve_y(
                blue if xv < 0.15 else coop_mask,
                x_axis,
                y_axis,
                float(xv),
                row_min,
                r_hi,
                pick_fn=lambda rs, c: min(rs, key=lambda s: s[0]),
                search_dx=22,
                note_prefix="Fig. 3H Coop. model (solid blue); ",
            )
            out.append(rd)
        return out

    if mode == "simple":
        # Dashed lower blue line: collect all dash pixel columns in bottom band [r_split, row_max]
        simple_mask = blue.copy()
        simple_mask[:r_split, :] = False
        # Compute per-column median row for every column containing blue dash pixels
        cols_ok, rows_ok, thicks = [], [], []
        for c in range(col_min, simple_mask.shape[1]):
            rs = C.runs_1d(simple_mask[r_split : row_max + 1, c])
            if rs:
                seg = max(rs, key=lambda t: t[1] - t[0])
                cols_ok.append(float(c) + 0.5)
                rows_ok.append(0.5 * (seg[0] + seg[1]) + r_split + 0.5)
                thicks.append(float(seg[1] - seg[0] + 1))
        cols_ok = np.array(cols_ok)
        rows_ok = np.array(rows_ok)
        med_thick = float(np.median(thicks)) if thicks else 12.0
        out = []
        for xv in at:
            c_f = float(x_axis.to_pixel(float(xv)))
            # Linear interpolation across nearest dash segments on left and right
            r_est = float(np.interp(c_f, cols_ok, rows_ok))
            half_px = max(0.5 * med_thick, 2.0) + 1.0 + y_axis.residual_px
            yv = float(y_axis.to_value(r_est))
            ya, yb = float(y_axis.to_value(r_est - half_px)), float(y_axis.to_value(r_est + half_px))
            out.append(
                {
                    "x": float(xv),
                    "y": yv,
                    "lo": min(ya, yb),
                    "hi": max(ya, yb),
                    "px": [c_f, r_est],
                    "px_lo": r_est - half_px,
                    "px_hi": r_est + half_px,
                    "note": f"Fig. 3H Simple model (dashed blue); cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
                }
            )
        return out
    raise ValueError(f"Unknown mode {mode!r}")


def read_series(image: np.ndarray, series: str, at: list[float]) -> list[dict]:
    if series == "data":
        return _read_data(image)
    if series in ("coop", "simple"):
        return _read_model_curve(image, series, at)
    raise ValueError(f"Unknown series {series!r}")
