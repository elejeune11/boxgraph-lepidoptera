"""Shared calibrated multi-curve reader for Main Figs. 8, 10A, and 12A (fast float32)."""
from __future__ import annotations

import numpy as np


def _clusters(idx: np.ndarray, gap: int = 8) -> list[np.ndarray]:
    if len(idx) == 0:
        return []
    splits = np.where(np.diff(idx) > gap)[0] + 1
    return list(np.split(idx, splits))


def _match_columns(
    sub_img: np.ndarray,
    sub_valid: np.ndarray,
    target_rgb: np.ndarray,
    other_rgbs: list[np.ndarray],
    rgb_tol: float,
    margin: float = 10.0,
) -> np.ndarray:
    """Fast vectorized swatch match on a (H, K, 3) slice."""
    arr = sub_img.astype(np.float32)
    diff = arr - target_rgb.astype(np.float32)
    d2_t = (diff * diff).sum(axis=2)
    mask = (d2_t <= rgb_tol * rgb_tol) & sub_valid
    if not mask.any():
        return mask
    d_t = np.sqrt(d2_t)
    for orgb in other_rgbs:
        od = arr - orgb.astype(np.float32)
        d_o = np.sqrt((od * od).sum(axis=2))
        mask &= (d_t < d_o - margin)
    return mask


def read_timeseries_curve(
    image: np.ndarray,
    at: list[float],
    swatch_idx: int,
    y_ticks_vals: list[float],
    x_ticks_vals: list[float] = (0.0, 500.0, 1000.0, 1500.0, 2000.0),
    rgb_tol: float = 42.0,
) -> list[dict]:
    """Calibrate axes from inward ticks, extract legend swatch RGBs from the image,
    and trace the requested curve across `at` with fast strided tracking."""
    H, W = image.shape[:2]
    dark = (image[:, :, 0] < 80) & (image[:, :, 1] < 80) & (image[:, :, 2] < 80)

    # 1. Plot box borders
    h_rows = np.where(dark.sum(axis=1) > 0.45 * W)[0]
    hr = [float(c.mean()) for c in _clusters(h_rows, gap=10)]
    r0, r1 = hr[0], hr[-1]
    v_cols = np.where(dark[int(r0):int(r1), :].sum(axis=0) > 0.85 * (r1 - r0))[0]
    vc = [float(c.mean()) for c in _clusters(v_cols, gap=10)]
    c0, c1 = vc[0], vc[-1]

    # 2. Inward ticks on left and bottom borders
    yt_mask = dark[int(r0 + 15):int(r1 - 15), int(c0 + 10):int(c0 + 28)].mean(axis=1) > 0.5
    yr = np.where(yt_mask)[0] + int(r0 + 15)
    y_int = [float(c.mean()) for c in _clusters(yr, gap=10)]
    y_tick_rows = np.array([r0] + y_int + [r1], dtype=float)
    y_vals_arr = np.array(y_ticks_vals, dtype=float)
    if len(y_tick_rows) != len(y_vals_arr):
        raise ValueError(f"Detected {len(y_tick_rows)} y-ticks, expected {len(y_vals_arr)}")
    yb, ya = np.polyfit(y_tick_rows, y_vals_arr, 1)
    y_res_px = float(np.max(np.abs(y_tick_rows - (y_vals_arr - ya) / yb)))

    xt_mask = dark[int(r1 - 28):int(r1 - 10), int(c0 + 15):int(c1 - 15)].mean(axis=0) > 0.5
    xc = np.where(xt_mask)[0] + int(c0 + 15)
    x_int = [float(c.mean()) for c in _clusters(xc, gap=10)]
    x_tick_cols = np.array([c0] + x_int + [c1], dtype=float)
    x_vals_arr = np.array(x_ticks_vals, dtype=float)
    if len(x_tick_cols) != len(x_vals_arr):
        raise ValueError(f"Detected {len(x_tick_cols)} x-ticks, expected {len(x_vals_arr)}")
    xb, xa = np.polyfit(x_tick_cols, x_vals_arr, 1)
    x_res_px = float(np.max(np.abs(x_tick_cols - (x_vals_arr - xa) / xb)))

    # 3. Locate legend box and read exact swatch RGBs from the image
    sr0, sr1 = int(r0 + 15), int(r0 + 0.60 * (r1 - r0))
    sc0, sc1 = int(c0 + 0.50 * (c1 - c0)), int(c1 - 10)
    leg_sub = dark[sr0:sr1, sc0:sc1]
    lh = np.where(leg_sub.sum(axis=1) > 0.15 * (c1 - c0))[0] + sr0
    lv = np.where(leg_sub.sum(axis=0) > 0.15 * (r1 - r0))[0] + sc0
    lr0, lr1 = int(lh.min()), int(lh.max())
    lc0, lc1 = int(lv.min()), int(lv.max())

    sw_c0, sw_c1 = lc0 + 30, lc0 + 150
    sw_patch = image[lr0 + 15:lr1 - 15, sw_c0:sw_c1, :]
    nonwhite = (sw_patch.mean(axis=2) < 230).mean(axis=1) > 0.5
    sw_rows = np.where(nonwhite)[0] + lr0 + 15
    sw_clusters = [c for c in _clusters(sw_rows, gap=10) if len(c) >= 4]
    swatches = []
    for cl in sw_clusters:
        rm = int(cl.mean())
        rgb = np.median(image[rm - 1:rm + 2, sw_c0 + 20:sw_c1 - 20, :].reshape(-1, 3), axis=0)
        swatches.append(rgb.astype(np.float32))
    if swatch_idx >= len(swatches):
        raise ValueError(f"Swatch {swatch_idx} out of range (found {len(swatches)})")
    target_rgb = swatches[swatch_idx]
    other_rgbs = [s for i, s in enumerate(swatches) if i != swatch_idx]
    is_black_target = (float(target_rgb.max()) < 40.0)
    is_green_target = (float(target_rgb[1]) > float(target_rgb[0]) + 30 and float(target_rgb[1]) > float(target_rgb[2]) + 30)

    # 4. Valid mask (excludes outer frame, legend box, and tick marks for black curve)
    ir0, ir1 = int(r0 + 8), int(r1 - 8)
    ic0, ic1 = int(c0 + 8), int(c1 - 8)
    valid = np.zeros((H, W), dtype=bool)
    valid[ir0:ir1, ic0:ic1] = True
    valid[max(0, lr0 - 8):min(H, lr1 + 8), max(0, lc0 - 8):min(W, lc1 + 8)] = False
    if is_black_target:
        for yt in y_tick_rows:
            valid[max(0, int(yt - 14)):min(H, int(yt + 14)), int(c0):int(c0 + 42)] = False
            valid[max(0, int(yt - 14)):min(H, int(yt + 14)), int(c1 - 42):int(c1)] = False
        for xt in x_tick_cols:
            valid[int(r0):int(r0 + 42), max(0, int(xt - 14)):min(W, int(xt + 14))] = False
            valid[int(r1 - 42):int(r1), max(0, int(xt - 14)):min(W, int(xt + 14))] = False

    # 5. Fast strided track extraction (every 8 px across the plot -> ~300-480 cols, <0.05s!)
    sample_cols = np.arange(ic0, ic1, 8)
    sub_img = image[ir0:ir1, sample_cols, :]
    sub_val = valid[ir0:ir1, sample_cols]
    sub_mask = _match_columns(sub_img, sub_val, target_rgb, other_rgbs, rgb_tol, margin=12.0)
    if is_green_target:
        # Guard against olive+black antialiasing edge pixels (which have R ≈ G)
        G_s = sub_img[:, :, 1].astype(np.float32)
        R_s = sub_img[:, :, 0].astype(np.float32)
        B_s = sub_img[:, :, 2].astype(np.float32)
        sub_mask &= (G_s - R_s > 25.0) & (G_s - B_s > 30.0)

    # Extract clean single-stroke track points from strided columns
    raw_c, raw_r, raw_h = [], [], []
    for k, col_px in enumerate(sample_cols):
        rr = np.where(sub_mask[:, k])[0]
        if len(rr) == 0:
            continue
        cls = [cl for cl in _clusters(rr, gap=14) if 3 <= len(cl) <= 450]
        if len(cls) == 1:
            cl = cls[0]
            raw_c.append(float(col_px))
            raw_r.append(float(cl.mean() + ir0))
            raw_h.append(float(len(cl)))
        elif len(cls) > 1 and raw_r:
            # Prefer cluster closest to recent track
            cl = min(cls, key=lambda cl_: abs(float(cl_.mean() + ir0) - raw_r[-1]))
            if abs(float(cl.mean() + ir0) - raw_r[-1]) < 140:
                raw_c.append(float(col_px))
                raw_r.append(float(cl.mean() + ir0))
                raw_h.append(float(len(cl)))

    raw_c = np.array(raw_c, dtype=float)
    raw_r = np.array(raw_r, dtype=float)
    raw_h = np.array(raw_h, dtype=float)

    # Remove isolated 1-point outlier spikes using median filter along raw_r
    if len(raw_r) >= 5:
        keep = np.ones(len(raw_r), dtype=bool)
        for j in range(len(raw_r)):
            lo_j, hi_j = max(0, j - 3), min(len(raw_r), j + 4)
            med = float(np.median(raw_r[lo_j:hi_j]))
            if abs(raw_r[j] - med) > 120.0:
                keep[j] = False
        raw_c, raw_r, raw_h = raw_c[keep], raw_r[keep], raw_h[keep]

    # 6. Evaluate at each requested x position (using exact narrow vertical slice!)
    out = []
    for xv in at:
        cx_f = (float(xv) - xa) / xb
        ci = int(round(cx_f))

        # Left-border (x = 0) extrapolation
        if abs(cx_f - c0) <= 12.0:
            near = (raw_c >= ic0) & (raw_c <= ic0 + 60)
            if near.any():
                idx_n = np.where(near)[0][:5]
                if len(idx_n) >= 2:
                    s_loc, i_loc = np.polyfit(raw_c[idx_n], raw_r[idx_n], 1)
                    cy = float(i_loc + s_loc * cx_f)
                else:
                    cy = float(raw_r[idx_n[0]])
                thick = float(np.median(raw_h[idx_n]))
                half_px = max(4.0, 0.5 * thick + 2.0 + y_res_px)
                y_val = float(ya + yb * cy)
                y_u = abs(yb) * half_px
                out.append({
                    "x": float(xv),
                    "y": round(y_val, 5),
                    "lo": round(y_val - y_u, 5),
                    "hi": round(y_val + y_u, 5),
                    "px": [round(cx_f, 2), round(cy, 2)],
                    "px_lo": round(cy - half_px, 2),
                    "px_hi": round(cy + half_px, 2),
                    "note": f"left-edge extrap; cal_res=({x_res_px:.2f},{y_res_px:.2f})px",
                })
                continue

        # Direct narrow-slice measurement at [ci-2 .. ci+2]
        sc_a = max(ic0, ci - 2)
        sc_b = min(ic1, ci + 3)
        found_direct = False
        if sc_a < sc_b:
            sl_img = image[ir0:ir1, sc_a:sc_b, :]
            sl_val = valid[ir0:ir1, sc_a:sc_b]
            sl_mask = _match_columns(sl_img, sl_val, target_rgb, other_rgbs, rgb_tol, margin=12.0)
            if is_green_target:
                G_s = sl_img[:, :, 1].astype(np.float32)
                R_s = sl_img[:, :, 0].astype(np.float32)
                B_s = sl_img[:, :, 2].astype(np.float32)
                sl_mask &= (G_s - R_s > 25.0) & (G_s - B_s > 30.0)
            rr = np.where(sl_mask.any(axis=1))[0]
            cls = [cl for cl in _clusters(rr, gap=14) if 3 <= len(cl) <= 450]
            if len(cls) >= 1:
                # Choose cluster consistent with nearby smooth track if multiple
                near_t = np.abs(raw_c - cx_f) <= 100
                if near_t.any():
                    ref_r = float(np.median(raw_r[near_t]))
                    cl = min(cls, key=lambda cl_: abs(float(cl_.mean() + ir0) - ref_r))
                    ok_track = abs(float(cl.mean() + ir0) - ref_r) <= 130.0
                else:
                    cl = max(cls, key=len)
                    ok_track = (len(cls) == 1)
                if ok_track:
                    cy = float(cl.mean() + ir0)
                    thick = float(len(cl))
                    half_px = max(3.0, 0.5 * thick + 1.5 + y_res_px)
                    y_val = float(ya + yb * cy)
                    y_u = abs(yb) * half_px
                    out.append({
                        "x": float(xv),
                        "y": round(y_val, 5),
                        "lo": round(y_val - y_u, 5),
                        "hi": round(y_val + y_u, 5),
                        "px": [round(cx_f, 2), round(cy, 2)],
                        "px_lo": round(cy - half_px, 2),
                        "px_hi": round(cy + half_px, 2),
                        "note": f"direct stroke h={thick:.0f}px; cal_res=({x_res_px:.2f},{y_res_px:.2f})px",
                    })
                    found_direct = True
        if found_direct:
            continue

        # Crossing interpolation or right-edge extrapolation from track
        left = np.where((raw_c < cx_f) & (raw_c >= cx_f - 200))[0]
        right = np.where((raw_c > cx_f) & (raw_c <= cx_f + 200))[0]
        if len(left) and len(right):
            il, ir_idx = left[-1], right[0]
            c_l, r_l = raw_c[il], raw_r[il]
            c_r, r_r = raw_c[ir_idx], raw_r[ir_idx]
            frac = (cx_f - c_l) / (c_r - c_l)
            cy = float(r_l + frac * (r_r - r_l))
            gap_px = c_r - c_l
            thick = float(0.5 * (raw_h[il] + raw_h[ir_idx]))
            half_px = max(5.0, 0.5 * thick + 0.15 * abs(r_r - r_l) + 0.05 * gap_px + y_res_px)
            y_val = float(ya + yb * cy)
            y_u = abs(yb) * half_px
            out.append({
                "x": float(xv),
                "y": round(y_val, 5),
                "lo": round(y_val - y_u, 5),
                "hi": round(y_val + y_u, 5),
                "px": [round(cx_f, 2), round(cy, 2)],
                "px_lo": round(cy - half_px, 2),
                "px_hi": round(cy + half_px, 2),
                "note": f"crossing interp gap={gap_px:.0f}px; cal_res=({x_res_px:.2f},{y_res_px:.2f})px",
            })
        elif len(left) and abs(cx_f - c1) <= 25:
            idx_l = left[-min(5, len(left)):]
            if len(idx_l) >= 2:
                s_loc, i_loc = np.polyfit(raw_c[idx_l], raw_r[idx_l], 1)
                cy = float(i_loc + s_loc * cx_f)
            else:
                cy = float(raw_r[idx_l[-1]])
            thick = float(np.median(raw_h[idx_l]))
            half_px = max(4.0, 0.5 * thick + 2.5 + y_res_px)
            y_val = float(ya + yb * cy)
            y_u = abs(yb) * half_px
            out.append({
                "x": float(xv),
                "y": round(y_val, 5),
                "lo": round(y_val - y_u, 5),
                "hi": round(y_val + y_u, 5),
                "px": [round(cx_f, 2), round(cy, 2)],
                "px_lo": round(cy - half_px, 2),
                "px_hi": round(cy + half_px, 2),
                "note": f"right-edge extrap; cal_res=({x_res_px:.2f},{y_res_px:.2f})px",
            })
        else:
            out.append({"x": float(xv), "y": None, "note": "occluded or not found"})

    return out
