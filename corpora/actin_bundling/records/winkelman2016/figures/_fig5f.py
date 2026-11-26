"""Calibrated reader for Fig. 5F (fascin domain length fraction vs [SNAP-647-Espin])."""

from __future__ import annotations

import numpy as np
from scipy import ndimage

import _calib as C


def _calibrate(image: np.ndarray):
    dark, gray, h_lines, v_lines = C.find_l_spines(
        image, thresh=100, row_frac=0.30, col_frac=0.30, row_range=(40, image.shape[0] - 80), col_range=(200, image.shape[1] - 20)
    )
    x_spine = max(h_lines, key=lambda t: t[0])
    y_spine = min(v_lines, key=lambda t: t[0])

    # Inward y ticks: 0.8, 0.6, 0.4, 0.2 (0.0 sits at x-spine center, 1.0 near top)
    yt_raw = C.find_ticks(dark, y_spine, (250, int(x_spine[1]) - 30), "y", side=+1, band=(4, 24), frac=0.35, max_width=45)
    yt = C.match_ticks_regular(yt_raw, 4, tol=8.0, min_step=250.0)
    y_axis = C.verify_axis("linear", list(zip(yt, [0.8, 0.6, 0.4, 0.2])), max_residual_px=2.0)

    # Inward x ticks: 0, 0.5, 1.0, 1.5, 2.0, 2.5 uM
    xt_raw = C.find_ticks(dark, x_spine, (int(y_spine[2]) + 15, image.shape[1] - 20), "x", side=-1, band=(4, 24), frac=0.35, max_width=45)
    xt = C.match_ticks_regular(xt_raw, 6, tol=8.0, min_step=230.0)
    x_axis = C.verify_axis("linear", list(zip(xt, [0.0, 0.5, 1.0, 1.5, 2.0, 2.5])), max_residual_px=3.0)
    return x_axis, y_axis, x_spine, y_spine


def read_series(image: np.ndarray, kind: str, at: list[float]) -> list[dict]:
    x_axis, y_axis, x_spine, y_spine = _calibrate(image)
    R = image[:, :, 0].astype(np.float32)
    G = image[:, :, 1].astype(np.float32)
    B = image[:, :, 2].astype(np.float32)

    if kind == "data":
        # Blue diamond markers (note: x=0 diamond overlaps y-spine!)
        blue = (B > 150) & (G > 90) & (B - R > 80)
        row_min = 20
        row_max = int(x_spine[1]) - 5
        col_min_blue = max(0, int(y_spine[1]) - 80)
        col_min = int(y_spine[2]) + 2
        blue[:row_min, :] = False
        blue[row_max + 1 :, :] = False
        blue[:, :col_min_blue] = False

        # Close across the black y-spine cutting the x=0 diamond so its center is exact!
        blue2 = blue[::2, ::2]
        blk2 = (R[::2, ::2] < 80) & (G[::2, ::2] < 80) & (B[::2, ::2] < 80)
        near_b2 = ndimage.binary_dilation(blue2, iterations=4)
        blue2_closed = ndimage.binary_closing(blue2 | (blk2 & near_b2), iterations=4)
        dist2 = ndimage.distance_transform_edt(blue2_closed)
        loc_max = (dist2 == ndimage.maximum_filter(dist2, size=15)) & (dist2 >= 12.0)
        lbl, n = ndimage.label(loc_max)
        pts = []
        H, W = blue.shape
        for i in range(1, n + 1):
            rr, cc = np.where(lbl == i)
            r0, c0 = float(rr.mean()) * 2.0, float(cc.mean()) * 2.0
            rad = int(round(dist2[int(round(rr.mean())), int(round(cc.mean()))] * 2.0))
            y0, y1 = max(0, int(r0) - rad), min(H, int(r0) + rad + 1)
            x0, x1 = max(0, int(c0) - rad), min(W, int(c0) + rad + 1)
            yy, xx = np.where(blue[y0:y1, x0:x1])
            d2 = (yy + y0 - r0) ** 2 + (xx + x0 - c0) ** 2
            keep = d2 <= (rad * 0.90) ** 2
            if keep.sum() > 20:
                pts.append((float((xx[keep] + x0).mean()) + 0.5, float((yy[keep] + y0).mean()) + 0.5))
            else:
                pts.append((c0 + 0.5, r0 + 0.5))
        pts.sort(key=lambda p: p[0])

        blk = (R < 80) & (G < 80) & (B < 80)
        blk[int(x_spine[1]) - 2 :, :] = False
        blk[:, :col_min] = False
        out = []
        for c_mid, r_mid in pts:
            xv = float(x_axis.to_value(c_mid))
            if abs(xv) < 0.025:
                xv = 0.0
            # Trace black vertical SEM error bar above and below blue diamond (stopping at first wide horizontal cap!)
            c0 = max(col_min, int(round(c_mid)) - 12)
            c1 = min(W, int(round(c_mid)) + 13)
            k0 = max(0, int(round(c_mid)) - 75)
            k1 = min(W, int(round(c_mid)) + 76)
            strip = blk[:, c0:c1].sum(axis=1) >= 4
            cap_w = blk[:, k0:k1].sum(axis=1)

            def trace_to_cap(start_r, stop_r, step):
                last_r = None
                gap = 0
                for r in range(start_r, stop_r, step):
                    if cap_w[r] >= 55:
                        # Found the horizontal cap! Center on contiguous cap rows
                        c_rows = [r]
                        rr = r + step
                        while 0 <= rr < H and cap_w[rr] >= 45 and abs(rr - r) <= 15:
                            c_rows.append(rr)
                            rr += step
                        return float(np.mean(c_rows)) + 0.5
                    if strip[r]:
                        last_r = r
                        gap = 0
                    else:
                        gap += 1
                        if gap > 35:
                            break
                return (float(last_r) + 0.5) if last_r is not None else None

            u_row = trace_to_cap(max(5, int(round(r_mid)) - 42), 4, -1)
            l_row = trace_to_cap(min(row_max, int(round(r_mid)) + 42), row_max + 1, +1)
            if u_row is not None and (r_mid - u_row) < 30:
                u_row = None
            if l_row is not None and (l_row - r_mid) < 30:
                l_row = None
            out.append(
                C.make_point_reading(
                    xv,
                    c_mid,
                    r_mid,
                    y_axis,
                    upper_err_row=u_row,
                    lower_err_row=l_row,
                    mirror_lower=True,
                    err_px=2.5 + y_axis.residual_px,
                    note=f"Fig. 5F data; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
                )
            )
        return out

    if kind == "fit":
        blk = (R < 80) & (G < 80) & (B < 80)
        row_min = 40
        row_max = int(x_spine[1]) - 10
        col_min = int(y_spine[2]) + 15
        blk[:row_min, :] = False
        blk[row_max + 1 :, :] = False
        blk[:, :col_min] = False
        return [
            C.sample_solid_curve_y(
                blk, x_axis, y_axis, float(xv), row_min, row_max, search_dx=12, note_prefix="Fig. 5F smooth fit; "
            )
            for xv in at
        ]
    raise ValueError(f"Unknown kind {kind!r}")
