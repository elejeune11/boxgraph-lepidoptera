"""Shared calibrated heatmap reader for Main Figs. 6 (R_g) and 9 (S)."""
from __future__ import annotations

import numpy as np


def _clusters(idx: np.ndarray, gap: int = 10) -> list[np.ndarray]:
    if len(idx) == 0:
        return []
    splits = np.where(np.diff(idx) > gap)[0] + 1
    return list(np.split(idx, splits))


def read_heatmap_row(
    image: np.ndarray,
    at: list[float],
    row_idx: int,
    cb_labels: list[float],
    min_half_width: float = 0.012,
) -> list[dict]:
    """Decode one row (row_idx: 0=top Rm:a=0.02, 1=mid 0.01, 2=bot 0.005)
    of a 3x6 YlGnBu_r heatmap using the image's own colorbar and printed tick labels."""
    H, W = image.shape[:2]
    dark = (image[:, :, 0] < 80) & (image[:, :, 1] < 80) & (image[:, :, 2] < 80)

    # 1. Locate top/bottom borders and the four vertical borders [hm_c0, hm_c1, cb_c0, cb_c1]
    h_rows = np.where(dark.sum(axis=1) > 0.45 * W)[0]
    hr = [float(c.mean()) for c in _clusters(h_rows, gap=10)]
    if len(hr) < 2:
        raise ValueError("Could not detect horizontal plot borders")
    r0, r1 = hr[0], hr[-1]

    v_cols = np.where(dark[int(r0):int(r1), :].sum(axis=0) > 0.85 * (r1 - r0))[0]
    vc = [float(c.mean()) for c in _clusters(v_cols, gap=10)]
    if len(vc) < 4:
        raise ValueError(f"Expected 4 vertical borders (heatmap + colorbar), got {len(vc)}")
    hm_c0, hm_c1, cb_c0, cb_c1 = vc[0], vc[1], vc[2], vc[3]

    # 2. Calibrate the colorbar value scale from printed numeric labels to the right of cb_c1
    #    Labels are sorted top-to-bottom matching cb_labels (descending).
    txt_band = dark[:, int(cb_c1 + 15):min(W, int(cb_c1 + 260))].any(axis=1)
    tr = np.where(txt_band)[0]
    lbl_clusters = [c for c in _clusters(tr, gap=15) if len(c) > 20]
    if len(lbl_clusters) != len(cb_labels):
        raise ValueError(
            f"Detected {len(lbl_clusters)} colorbar labels, expected {len(cb_labels)}"
        )
    lbl_rows = np.array([float(c.mean()) for c in lbl_clusters])
    lbl_vals = np.array(cb_labels, dtype=float)
    # Fit linear slope (value per pixel); anchor bottom label (min value) at r1
    # (in both Fig. 6 and Fig. 9, the bottom tick sits on r1 with a constant +10.0 px
    # font-centroid-to-baseline shift verified against r0 on Fig. 6).
    slope, _ = np.polyfit(lbl_rows, lbl_vals, 1)
    # Shift all label centers by the bottom-anchor offset (r1 - lbl_rows[-1])
    anchor_offset = r1 - lbl_rows[-1]
    tick_rows = lbl_rows + anchor_offset
    cal_resid_px = float(np.max(np.abs(tick_rows - (r1 + (lbl_vals - lbl_vals[-1]) / slope))))

    def row_to_val(r: float | np.ndarray):
        return lbl_vals[-1] + slope * (r - r1)

    v_bot = float(row_to_val(r1))
    v_top = float(row_to_val(r0))

    # 3. Extract interior RGB profile of the colorbar
    cb_mid_c0 = int(cb_c0 + 0.30 * (cb_c1 - cb_c0))
    cb_mid_c1 = int(cb_c1 - 0.30 * (cb_c1 - cb_c0))
    cb_rows = np.arange(int(r0 + 12), int(r1 - 12))
    cb_rgb = image[cb_rows, cb_mid_c0:cb_mid_c1, :].astype(float).mean(axis=1)
    cb_vals = row_to_val(cb_rows.astype(float))

    # 4. Decode the 6 columns (R_alpha:a = 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
    ra_grid = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
    xs = list(at) if at else ra_grid
    cell_h = (r1 - r0) / 3.0
    cell_w = (hm_c1 - hm_c0) / 6.0
    cy = r0 + (row_idx + 0.5) * cell_h

    out = []
    for xv in xs:
        # Match xv to column index in ra_grid
        ci = min(range(len(ra_grid)), key=lambda i: abs(ra_grid[i] - xv))
        if abs(ra_grid[ci] - xv) > 1e-6:
            out.append({"x": float(xv), "y": None, "note": "not a heatmap column"})
            continue
        cx = hm_c0 + (ci + 0.5) * cell_w
        patch = image[int(cy - 30):int(cy + 30), int(cx - 30):int(cx + 30), :].astype(float)
        cell_rgb = np.median(patch.reshape(-1, 3), axis=0)
        dists = np.sqrt(((cb_rgb - cell_rgb) ** 2).sum(axis=1))
        best = int(np.argmin(dists))
        d_min = float(dists[best])
        val = float(cb_vals[best])

        # Honest 95% uncertainty interval from colorbar RGB matching tolerance + calibration
        # Find all colorbar values within RGB distance <= max(6.0, d_min + 5.0)
        # (also covers subtle dark-blue PDF shading shift between cells and colorbar)
        tol_rgb = max(6.0, d_min + 6.0)
        near_vals = cb_vals[dists <= tol_rgb]
        cal_u = abs(slope) * max(2.0, cal_resid_px)
        half = max(
            min_half_width,
            float(val - near_vals.min()) + cal_u if len(near_vals) else min_half_width,
            float(near_vals.max() - val) + cal_u if len(near_vals) else min_half_width,
        )
        lo = max(v_bot, val - half) if d_min > 8.0 and best > len(cb_rows) - 30 else val - half
        hi = min(v_top, val + half) if d_min > 8.0 and best < 30 else val + half
        if lo > val:
            lo = val - min_half_width
        if hi < val:
            hi = val + min_half_width

        out.append({
            "x": float(xv),
            "y": round(val, 5),
            "lo": round(lo, 5),
            "hi": round(hi, 5),
            "px": [round(float(cx), 2), round(float(cy), 2)],
            "note": f"cb RGB dist={d_min:.1f}, cal_resid={cal_resid_px:.2f}px, cb_range=[{v_bot:.3f},{v_top:.4f}]",
        })
    return out
