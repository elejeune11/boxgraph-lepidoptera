"""Calibrated reader for SI Fig. S1 (PDF p. 18 / SI-9): normalized relative actin density rho/rho_max
for the 5 experimental samples (3 'low' samples indexed x=1,2,3 and 2 'high' samples indexed x=4,5)."""
from __future__ import annotations

import numpy as np
from scipy import ndimage
from boxgraph import digitize as D


def _clusters(mask: np.ndarray, min_len: int = 1, gap: int = 2) -> list[tuple[int, int, float]]:
    idx = np.where(mask)[0]
    if len(idx) == 0:
        return []
    splits = np.where(np.diff(idx) > gap)[0] + 1
    out = []
    for g in np.split(idx, splits):
        if len(g) >= min_len:
            out.append((int(g[0]), int(g[-1]), float(g.mean())))
    return out


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    gray = image.astype(float).mean(axis=2)
    black = gray < 130.0
    H, W = gray.shape

    # Locate plot spines
    vsp = _clusters(black.mean(axis=0) > 0.55)
    hsp = _clusters(black.mean(axis=1) > 0.55)
    c0, c1 = vsp[0][2], vsp[-1][2]
    r0, r1 = hsp[0][2], hsp[-1][2]
    c0i, c1i, r0i, r1i = int(round(c0)), int(round(c1)), int(round(r0)), int(round(r1))

    # Bottom inward ticks: "low" and "high" columns
    xt_cls = _clusters(black[r1i - 45:r1i - 15, c0i + 30:c1i - 30].mean(axis=0) > 0.70)
    if len(xt_cls) != 2:
        raise RuntimeError(f"Fig. S1: expected 2 x ticks (low, high), got {xt_cls}")
    c_low = xt_cls[0][2] + c0i + 30
    c_high = xt_cls[1][2] + c0i + 30

    # Left inward y ticks: 1.0, 0.8, 0.6, 0.4, 0.2 (excluding horizontal dashed line at 0.5 which spans full width!)
    yt_cls = _clusters(black[r0i + 30:r1i - 30, c0i + 15:c0i + 45].mean(axis=1) > 0.70)
    # Exclude the dashed line (which also has black pixels near the right side of the plot!)
    yt = []
    for a, b, m in yt_cls:
        rm = int(round(m + r0i + 30))
        # A genuine tick has white paper to its right (cols c0+120..c0+350), whereas dashed line has ink there
        if black[rm - 3:rm + 4, c0i + 120:c0i + 350].mean() < 0.08:
            yt.append(m + r0i + 30)
    if len(yt) != 5:
        raise RuntimeError(f"Fig. S1: expected 5 y ticks (1.0, 0.8, 0.6, 0.4, 0.2), got {yt}")
    ya = D.Axis("linear", [(yt[0], 1.0), (yt[1], 0.8), (yt[2], 0.6), (yt[3], 0.4), (yt[4], 0.2)])

    # Detect the 5 open circles inside the plot interior using enclosed white-interior holes!
    # Mask out the horizontal dashed cutoff line outside the circles, or simply detect circle rings in each column:
    # Each open circle has outer diameter ~124 px (radius ~62 px) and thick black ring (~14 px).
    # In the 'low' and 'high' vertical strips (c_col +/- 90 px), remove thin horizontal dashed line (thickness <= 12 px
    # spanning full strip width) before finding white enclosed interiors!
    interior = np.zeros((r1i - r0i - 30, c1i - c0i - 30), dtype=bool)
    roi_black = black[r0i + 15:r1i - 15, c0i + 15:c1i - 15].copy()
    # Erase dashed line pixels at rows where > 35% of the full plot width is black (purely the 0.5 cutoff line)
    dash_rows = np.where(roi_black.mean(axis=1) > 0.30)[0]
    # Only erase dashed line OUTSIDE the circle ring (actually inside the low circle it bridges horizontally;
    # erasing dashed rows splits that ring by ~8 px vertically, which a 10-px vertical closing restores cleanly!)
    roi_no_dash = roi_black.copy()
    for dr in dash_rows:
        # Erase horizontal segments NOT part of the circle ring: easiest is to fill enclosed white holes on
        # roi_black directly (for circle at 0.5, its white interior splits into top half-moon & bottom half-moon,
        # which merge into one circle center if we take connected components of (filled - outer_bg) within 140 px!).
        pass

    # Flood-fill outer white background of roi_black: anything white NOT reachable from border is inside a circle!
    white = ~roi_black
    # Label all white regions
    lbl, nlbl = ndimage.label(white)
    # Border labels are background
    bg_labels = set(np.unique(lbl[0, :])) | set(np.unique(lbl[-1, :])) | set(np.unique(lbl[:, 0])) | set(np.unique(lbl[:, -1]))
    holes = np.zeros_like(white)
    for k in range(1, nlbl + 1):
        if k in bg_labels:
            continue
        area = int((lbl == k).sum())
        if 300 < area < 20000:
            holes |= (lbl == k)

    # Dilate holes vertically by 10 px ONLY to merge the two half-moons of the circle bisected by the 0.5 dashed line
    # (note: the two touching circles at the bottom are separated by ~28 px of black ring, so 10-px dilation keeps them separate!)
    # Check: 28 px black junction -> half is 14 px, so a 7-px vertical dilation bridges the 8-px dashed line while keeping the 28-px junction separate!
    struct_v = np.ones((11, 1), dtype=bool)
    merged = ndimage.binary_dilation(holes, structure=struct_v)
    mlbl, ncomp = ndimage.label(merged)

    low_pts, high_pts = [], []
    for k in range(1, ncomp + 1):
        # Compute centroid on original un-dilated hole pixels belonging to component k
        rr, cc = np.nonzero(holes & (mlbl == k))
        if len(rr) < 500:
            continue
        r_cent = float(rr.mean()) + r0i + 15
        c_cent = float(cc.mean()) + c0i + 15
        val = ya.to_value(r_cent)
        if abs(c_cent - c_low) < abs(c_cent - c_high):
            low_pts.append((val, c_cent, r_cent, "low"))
        else:
            high_pts.append((val, c_cent, r_cent, "high"))

    low_pts.sort(key=lambda t: t[0])
    high_pts.sort(key=lambda t: t[0])
    all_pts = low_pts + high_pts
    if len(all_pts) != 5:
        raise RuntimeError(f"Fig. S1: expected 5 sample circles, found {len(all_pts)}: {all_pts}")

    out = []
    half_px = max(4.0, ya.residual_px + 3.0)
    for idx, (val, c_cent, r_cent, grp) in enumerate(all_pts, start=1):
        hi = ya.to_value(r_cent - half_px)
        lo = ya.to_value(r_cent + half_px)
        out.append({
            "x": float(idx),
            "y": float(val),
            "lo": float(lo),
            "hi": float(hi),
            "px": [float(c_cent), float(r_cent)],
            "px_lo": float(r_cent - half_px),
            "px_hi": float(r_cent + half_px),
            "note": f"group={grp} y_res={ya.residual_px:.2f}px",
        })
    return out
