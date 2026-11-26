"""Shared calibrated reader for S4 Text Excel vector plots (S4 Figs. 1-5), fast version."""
from __future__ import annotations

import numpy as np
from scipy import ndimage


def _clusters(idx: np.ndarray, gap: int = 8) -> list[np.ndarray]:
    if len(idx) == 0:
        return []
    splits = np.where(np.diff(idx) > gap)[0] + 1
    return list(np.split(idx, splits))


def _series_mask(sub: np.ndarray, color: str) -> np.ndarray:
    R = sub[:, :, 0].astype(np.int16)
    G = sub[:, :, 1].astype(np.int16)
    B = sub[:, :, 2].astype(np.int16)
    if color == "blue":
        return (B - R > 30) & (B > 120)
    if color == "red":
        return (R - B > 35) & (R - G > 25) & (R > 120)
    if color == "green":
        return (G - R > 12) & (G - B > 35) & (G > 120)
    raise ValueError(f"Unknown Excel series color {color!r}")


def _disk(radius: int) -> np.ndarray:
    yy, xx = np.ogrid[-radius:radius + 1, -radius:radius + 1]
    return (xx * xx + yy * yy) <= radius * radius


_DISK_8 = _disk(8)
_DISK_6 = _disk(6)


def _find_marker_center(mask: np.ndarray, rel_ci: int, sr0: int, sc0: int) -> tuple[float, float, float] | None:
    """Use morphological opening to isolate the marker body from steep connecting lines and thin caps."""
    H_s, W_s = mask.shape
    c0, c1 = max(0, rel_ci - 26), min(W_s, rel_ci + 27)
    sub = mask[:, c0:c1]
    for struct in (_DISK_8, _DISK_6):
        opened = ndimage.binary_opening(sub, structure=struct)
        lbl, n_comp = ndimage.label(opened)
        if n_comp > 0:
            best = None
            best_score = -1.0
            for cid in range(1, n_comp + 1):
                rr, cc = np.nonzero(lbl == cid)
                mc = float(cc.mean()) + c0
                mr = float(rr.mean())
                dist_x = abs(mc - rel_ci)
                if dist_x > 18.0:
                    continue
                score = float(len(rr)) / (1.0 + 0.1 * dist_x)
                if score > best_score:
                    best_score = score
                    best = (mc + sc0, mr + sr0, float(rr.max() - rr.min() + 1))
            if best is not None:
                return best
    return None


def read_s4_excel(
    image: np.ndarray,
    at: list[float],
    series: str,
    y_grid_vals: list[float],
    x_grid_vals: list[float],
    legend_frac_box: tuple[float, float, float, float],
    has_error_bars: bool = False,
    fallback_on_red_xs: tuple[float, ...] = (),
    fallback_on_green_xs: tuple[float, ...] = (),
    bridge_other_color: str | None = None,
    non_negative: bool = False,
) -> list[dict]:
    """Calibrate an S4 Text Excel chart from its horizontal gridlines and x-axis labels,
    then extract exact marker centers (and symmetric error-bar caps if present)."""
    H, W = image.shape[:2]

    # 1. Fast horizontal gridline detection using strided columns across [0.15W..0.95W]
    sub_g = image[:, int(0.15 * W):int(0.95 * W):16, :]
    g_max = sub_g.max(axis=2)
    g_min = sub_g.min(axis=2)
    g_mean = (sub_g[:, :, 0].astype(np.int16) + sub_g[:, :, 1].astype(np.int16) + sub_g[:, :, 2].astype(np.int16)) // 3
    grid_s = (g_mean < 215) & ((g_max - g_min) < 15)
    h_rows = np.where(grid_s.mean(axis=1) > 0.72)[0]
    hcl = [float(c.mean()) for c in _clusters(h_rows, gap=5)]
    grid_rows = np.array([r for r in hcl if 280.0 < r < 3400.0], dtype=float)
    y_vals_arr = np.array(y_grid_vals, dtype=float)
    if len(grid_rows) != len(y_vals_arr):
        raise ValueError(f"Detected {len(grid_rows)} horizontal gridlines, expected {len(y_vals_arr)}")
    yb, ya = np.polyfit(grid_rows, y_vals_arr, 1)
    y_res_px = float(np.max(np.abs(grid_rows - (y_vals_arr - ya) / yb)))
    r_top, r_bot = float(grid_rows[0]), float(grid_rows[-1])

    # 2. Horizontal span of gridlines from top gridline (5-row band)
    rt_i = int(round(r_top))
    top_band = image[rt_i - 2:rt_i + 3, :, :]
    tb_mean = (top_band[:, :, 0].astype(np.int16) + top_band[:, :, 1].astype(np.int16) + top_band[:, :, 2].astype(np.int16)) // 3
    tb_grid = (tb_mean < 215) & ((top_band.max(axis=2) - top_band.min(axis=2)) < 15)
    gcols = np.where(tb_grid.any(axis=0))[0]
    main_seg = max(_clusters(gcols, gap=10), key=len)
    gx0, gx1 = float(main_seg[0]), float(main_seg[-1])

    # 3. X-axis category centers from tick labels directly below r_bot inside [gx0+10, gx1-10]
    rb_i = int(round(r_bot))
    lbl_band = image[rb_i + 35:rb_i + 125, int(gx0 + 10):int(gx1 - 10), :]
    lbl_gray = (lbl_band[:, :, 0].astype(np.int16) + lbl_band[:, :, 1].astype(np.int16) + lbl_band[:, :, 2].astype(np.int16)) // 3
    xl_cols = np.where((lbl_gray < 150).any(axis=0))[0] + int(gx0 + 10)
    x_cl = [c for c in _clusters(xl_cols, gap=45) if len(c) >= 10]
    if len(x_cl) != len(x_grid_vals):
        n_cat = len(x_grid_vals)
        x_centers = np.array([gx0 + (i + 0.5) * (gx1 - gx0) / n_cat for i in range(n_cat)], dtype=float)
        x_res_px = 2.0
    else:
        raw_xc = np.array([float(c.mean()) for c in x_cl], dtype=float)
        idx_arr = np.arange(len(x_grid_vals), dtype=float)
        xb_i, xa_i = np.polyfit(idx_arr, raw_xc, 1)
        x_centers = xa_i + xb_i * idx_arr
        x_res_px = float(np.max(np.abs(raw_xc - x_centers)))

    # 4. Legend rectangle in plot coordinates
    lf_c0, lf_r0, lf_c1, lf_r1 = legend_frac_box
    leg_c0 = int(gx0 + lf_c0 * (gx1 - gx0))
    leg_c1 = int(gx0 + lf_c1 * (gx1 - gx0))
    leg_r0 = int(r_top + lf_r0 * (r_bot - r_top))
    leg_r1 = int(r_top + lf_r1 * (r_bot - r_top))

    sr0 = max(0, int(r_top - 40))
    sr1 = min(H, int(r_bot + 50))
    xs = list(at) if at else list(x_grid_vals)

    out = []
    for xv in xs:
        ci_idx = min(range(len(x_grid_vals)), key=lambda i: abs(x_grid_vals[i] - xv))
        if abs(x_grid_vals[ci_idx] - xv) > 1e-6:
            out.append({"x": float(xv), "y": None, "note": "not an x category"})
            continue
        cx_exp = float(x_centers[ci_idx])
        ci = int(round(cx_exp))

        w_half = 44
        sc0 = max(int(gx0 - 10), ci - w_half)
        sc1 = min(int(gx1 + 10), ci + w_half + 1)
        strip_img = image[sr0:sr1, sc0:sc1, :]
        valid_s = np.ones((sr1 - sr0, sc1 - sc0), dtype=bool)
        lr_a, lr_b = max(0, leg_r0 - sr0), min(sr1 - sr0, leg_r1 - sr0)
        lc_a, lc_b = max(0, leg_c0 - sc0), min(sc1 - sc0, leg_c1 - sc0)
        if lr_a < lr_b and lc_a < lc_b:
            valid_s[lr_a:lr_b, lc_a:lc_b] = False

        mask_t = _series_mask(strip_img, series) & valid_s
        rel_ci = ci - sc0
        found_m = _find_marker_center(mask_t, rel_ci, sr0, sc0)

        # Fallback if completely covered by a coincident top-drawn series (red or green)
        occluded_by = None
        if found_m is None and any(abs(xv - fx) < 1e-6 for fx in fallback_on_red_xs):
            fb_m = _series_mask(strip_img, "red") & valid_s
            found_fb = _find_marker_center(fb_m, rel_ci, sr0, sc0)
            if found_fb is not None:
                occluded_by = "red"
                cx, cy, _ = found_fb
                half_px = 12.0 + y_res_px
                note_extra = "coincident/occluded behind red marker"
                active_m = fb_m
        elif found_m is None and any(abs(xv - fx) < 1e-6 for fx in fallback_on_green_xs):
            fb_m = _series_mask(strip_img, "green") & valid_s
            found_fb = _find_marker_center(fb_m, rel_ci, sr0, sc0)
            if found_fb is not None:
                occluded_by = "green"
                cx, cy, _ = found_fb
                half_px = 12.0 + y_res_px
                note_extra = "coincident/occluded behind green marker"
                active_m = fb_m

        if occluded_by is None:
            if found_m is not None:
                cx, cy, mh = found_m
                half_px = max(3.5, y_res_px + 2.5)
                note_extra = f"marker h={mh:.0f}px"
                active_m = mask_t
            else:
                # Partially covered marker: interpolate connecting line on left & right
                l_rr = np.where(mask_t[:, max(0, rel_ci - 36):max(0, rel_ci - 18)].any(axis=1))[0]
                r_rr = np.where(mask_t[:, min(sc1 - sc0, rel_ci + 18):min(sc1 - sc0, rel_ci + 37)].any(axis=1))[0]
                if len(l_rr) and len(r_rr):
                    cy = 0.5 * (float(np.median(l_rr)) + float(np.median(r_rr))) + sr0
                    half_px = 6.0 + y_res_px
                    note_extra = "interpolated from adjacent line segments"
                elif len(l_rr):
                    cy = float(np.median(l_rr)) + sr0
                    half_px = 8.0 + y_res_px
                    note_extra = "from left adjacent line segment"
                else:
                    rr_any = np.where(mask_t[:, max(0, rel_ci - 10):min(sc1 - sc0, rel_ci + 11)].any(axis=1))[0]
                    if len(rr_any) == 0:
                        out.append({"x": float(xv), "y": None, "note": "marker not found"})
                        continue
                    cy = float(np.median(rr_any)) + sr0
                    half_px = 10.0 + y_res_px
                    note_extra = "partial visible marker"
                cx = cx_exp
                active_m = mask_t

        # Error bars (if present on this series)
        bar_fields = {}
        px_lo, px_hi = cy - half_px, cy + half_px
        if has_error_bars:
            cap_c0 = max(0, rel_ci - 40)
            cap_c1 = min(sc1 - sc0, rel_ci + 41)
            eb_m = mask_t if mask_t[:, cap_c0:cap_c1].sum() >= 60 else active_m
            eb_row_w = eb_m[:, cap_c0:cap_c1].sum(axis=1)
            stem_self = eb_m[:, max(0, rel_ci - 6):min(sc1 - sc0, rel_ci + 7)].any(axis=1)
            if bridge_other_color is not None:
                other_m = _series_mask(strip_img, bridge_other_color) & valid_s
                stem_bridge = stem_self | other_m[:, max(0, rel_ci - 6):min(sc1 - sc0, rel_ci + 7)].any(axis=1)
            else:
                stem_bridge = stem_self
            rel_cy = int(round(cy - sr0))

            # Walk up and down along stem, tracking furthest self-colored stem pixel
            r_up_self = rel_cy
            gap = 0
            for r_i in range(rel_cy, -1, -1):
                if stem_bridge[r_i] or abs(r_i - rel_cy) <= 16:
                    if stem_self[r_i]:
                        r_up_self = r_i
                    gap = 0
                else:
                    gap += 1
                    if gap > 8:
                        break
            r_dn_self = rel_cy
            gap = 0
            for r_i in range(rel_cy, len(stem_bridge)):
                if stem_bridge[r_i] or abs(r_i - rel_cy) <= 16:
                    if stem_self[r_i]:
                        r_dn_self = r_i
                    gap = 0
                else:
                    gap += 1
                    if gap > 8:
                        break

            # Detect TRUE horizontal error-bar caps: wide (>= 28 px) AND thin (<= 11 px tall)
            # AND strictly outside the marker radius (>= 16 px from rel_cy)!
            marker_r = 14.0
            cap_rows = np.where(eb_row_w >= 28)[0]
            cap_cls = [c for c in _clusters(cap_rows, gap=3) if len(c) <= 12]
            top_caps = [float(c.mean()) for c in cap_cls if (r_up_self - 6) <= float(c.mean()) <= (rel_cy - marker_r - 2.0)]
            bot_caps = [float(c.mean()) for c in cap_cls if (rel_cy + marker_r + 2.0) <= float(c.mean()) <= (r_dn_self + 6)]

            has_top_cap = len(top_caps) > 0
            has_bot_cap = len(bot_caps) > 0
            top_px = (min(top_caps) + sr0) if has_top_cap else float(r_up_self + sr0)
            bot_px = (max(bot_caps) + sr0) if has_bot_cap else float(r_dn_self + sr0)
            up_len = cy - top_px
            dn_len = bot_px - cy

            if not has_top_cap and not has_bot_cap:
                top_px = cy - marker_r
                bot_px = cy + marker_r
                note_extra += "; eb hidden inside marker (<= marker radius)"
            elif has_top_cap and (not has_bot_cap or (bridge_other_color and dn_len < 0.70 * up_len)):
                bot_px = cy + up_len
                note_extra += f"; lower cap obscured by {bridge_other_color or 'overlap'}, mirrored from upper cap"
            elif has_bot_cap and (not has_top_cap or (bridge_other_color and up_len < 0.70 * dn_len)):
                top_px = cy - dn_len
                note_extra += f"; upper cap obscured by {bridge_other_color or 'overlap'}, mirrored from lower cap"
            else:
                note_extra += f"; eb caps=[{top_px:.1f},{bot_px:.1f}]"

            bar_hi = float(ya + yb * top_px)
            bar_lo = float(ya + yb * bot_px)
            bar_lo, bar_hi = min(bar_lo, bar_hi), max(bar_lo, bar_hi)
            if non_negative:
                bar_lo = max(0.0, bar_lo)
            bar_u = abs(yb) * max(3.5, y_res_px + 2.5)
            bar_lo_outer = max(0.0, bar_lo - bar_u) if non_negative else (bar_lo - bar_u)
            bar_fields = {
                "bar_lo": round(bar_lo, 6),
                "bar_hi": round(bar_hi, 6),
                "bar_lo_outer": round(bar_lo_outer, 6),
                "bar_hi_outer": round(bar_hi + bar_u, 6),
            }
            px_lo, px_hi = top_px, bot_px

        y_val = float(ya + yb * cy)
        y_u = abs(yb) * half_px
        lo_v = max(0.0, y_val - y_u) if (non_negative and y_val >= 0.0) else (y_val - y_u)
        if y_val < 0.0 and non_negative:
            y_val = 0.0
            lo_v = 0.0
        rec = {
            "x": float(xv),
            "y": round(y_val, 6),
            "lo": round(lo_v, 6),
            "hi": round(y_val + y_u, 6),
            "px": [round(cx, 2), round(cy, 2)],
            "px_lo": round(px_lo, 2),
            "px_hi": round(px_hi, 2),
            "note": f"{note_extra}; cal_res=({x_res_px:.2f},{y_res_px:.2f})px",
        }
        rec.update(bar_fields)
        out.append(rec)

    return out
