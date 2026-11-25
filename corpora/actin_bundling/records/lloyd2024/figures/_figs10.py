"""Shared calibrated reader for SI Fig. 10A (wavelength spectra) and 10B (angular intensity) 6x6 facet grids (PDF p. 34, 900 dpi, rotated 90 deg CW)."""

from __future__ import annotations

import numpy as np
from boxgraph import digitize as D
from _common import find_runs, run_centers, y_interval


def read_figs10_peaks(image: np.ndarray, panel: str, group: str) -> list[dict]:
    """Rotate SI Fig. 10A or 10B 90 deg CCW upright and extract peak reflectance (%) per individual facet.

    group == 'control': first 15 facets (all red GIM_TREAT_* individuals, x=1..15 in row-major facet order).
    group == 'treated': remaining 21 facets (all teal TREATED_* individuals, x=1..21 in row-major facet order).
    """
    is_dots = (panel == "10B")
    W_orig = image.shape[1]
    rot = np.rot90(image, k=1)

    R = rot[:, :, 0].astype(float)
    G = rot[:, :, 1].astype(float)
    B = rot[:, :, 2].astype(float)
    gray = (R + G + B) / 3.0

    # Locate the 6 grey facet header rows and 6 facet column intervals
    strip = (np.abs(R - G) < 8.0) & (np.abs(G - B) < 8.0) & (gray > 190.0) & (gray < 230.0)
    hdrs = [(s, e) for s, e, _ in find_runs(np.where(strip.mean(axis=1) > 0.25)[0], gap=10) if (e - s) > 80]
    cols = [(s, e) for s, e, _ in find_runs(np.where(strip.mean(axis=0) > 0.08)[0], gap=10) if (e - s) > 400]
    if len(hdrs) != 6 or len(cols) != 6:
        raise RuntimeError(f"SI Fig. {panel}: expected 6x6 facets, found {len(hdrs)}x{len(cols)}")

    c0_left = cols[0][0]
    red_m = (R - G > 50.0) & (R - B > 50.0) & (R > 150.0)
    teal_m = (B - R > 40.0) & (G - R > 40.0) & (G > 130.0)
    y_lbls = [25.0, 20.0, 15.0, 10.0, 5.0, 0.0] if is_dots else [25.0, 20.0, 15.0, 10.0, 5.0]

    out = []
    idx_out = 0
    for ri, (hs, he) in enumerate(hdrs):
        rb_end = hdrs[ri + 1][0] - 25 if ri + 1 < len(hdrs) else rot.shape[0] - 180
        # Top tick (25%) sits flush with the header bottom edge he; search from he - 6 downward
        t_rows = run_centers(np.where(gray[he - 6:rb_end, c0_left - 10] < 160.0)[0] + he - 6, gap=4)
        if len(t_rows) != len(y_lbls):
            raise RuntimeError(f"SI Fig. {panel} row {ri}: expected {len(y_lbls)} y-ticks, found {t_rows}")
        ya = D.Axis("linear", list(zip(t_rows, y_lbls)))

        for ci, (cs, ce) in enumerate(cols):
            sub_r = red_m[he + 8:rb_end - 12, cs + 10:ce - 10]
            sub_t = teal_m[he + 8:rb_end - 12, cs + 10:ce - 10]
            nr, nt = int(sub_r.sum()), int(sub_t.sum())
            facet_group = "control" if nr > nt else "treated"
            if facet_group != group:
                continue
            idx_out += 1
            m_use = sub_r if facet_group == "control" else sub_t
            if not is_dots:
                # Smooth spectrum line: vertical centerline per column, pick highest (minimum row) peak
                c_ok = np.where(m_use.any(axis=0))[0]
                r_cents = np.array([np.where(m_use[:, c])[0].mean() + he + 8 + 0.5 for c in c_ok])
                k_best = int(np.argmin(r_cents))
                r_pk = float(r_cents[k_best])
                c_pk = float(c_ok[k_best] + cs + 10 + 0.5)
                half_px = 4.0
            else:
                # 11 discrete circular markers: centroid of highest marker (smallest row)
                rr, cc = np.where(m_use)
                c_top = int(cc[int(np.argmin(rr))])
                near = (np.abs(cc - c_top) <= 24) & (rr <= int(rr.min()) + 50)
                r_pk = float(rr[near].mean() + he + 8 + 0.5)
                c_pk = float(cc[near].mean() + cs + 10 + 0.5)
                half_px = 4.0

            y_val, lo, hi, _, _ = y_interval(ya, r_pk, half_px)
            # Inverse of np.rot90(img, k=1): rot[r_rot, c_rot] == img[c_rot, W_orig - 1 - r_rot]
            col_orig = float(W_orig - 1.0 - r_pk)
            row_orig = float(c_pk)
            out.append({
                "x": float(idx_out),
                "y": float(y_val), "lo": float(lo), "hi": float(hi),
                "px": [col_orig, row_orig],
                "note": f"SI Fig. {panel} {group} #{idx_out} (grid r{ri+1}c{ci+1}, red={nr}, teal={nt}, ry={ya.residual_px:.2f} px)",
            })

    expected = 15 if group == "control" else 21
    if len(out) != expected:
        raise RuntimeError(f"SI Fig. {panel} ({group}): expected {expected} facets, got {len(out)}")
    return out
