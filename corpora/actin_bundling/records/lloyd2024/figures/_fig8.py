"""Shared readers for Fig. 8E reflectance spectra and Fig. 8N-O dot+mean/SD plots (PDF p. 11, 1200 dpi)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage
from boxgraph import digitize as D
from _common import find_runs, run_centers, y_interval


def read_fig8e_spectrum(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    """Read Control (red) or Cyto D (blue) mean reflectance curve and +-SE shaded ribbon from Fig. 8E."""
    gray = image.astype(float).mean(axis=2)
    dark = gray < 120.0
    H, W = dark.shape

    c_lines = run_centers(np.where(dark[120:1940, :].mean(axis=0) > 0.70)[0], gap=4)
    c0 = [c for c in c_lines if c > 200.0][0]
    r_lines = run_centers(np.where(dark[:, int(c0) + 30:].mean(axis=1) > 0.70)[0], gap=4)
    r0, r1 = r_lines[-2], r_lines[-1]

    yt = run_centers(np.where(dark[int(r0) + 10:int(r1) - 10, int(round(c0)) - 10])[0] + int(r0) + 10, gap=4)
    xt = run_centers(np.where(dark[int(round(r1)) + 10, int(round(c0)) + 10:W - 10])[0] + int(round(c0)) + 10, gap=4)
    if len(xt) != 8 or len(yt) != 5:
        raise RuntimeError(f"Fig. 8E expected 8 x-ticks and 5 y-ticks, got xt={xt}, yt={yt}")

    xa = D.Axis("linear", list(zip(xt, [350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0])))
    ya = D.Axis("linear", list(zip(yt, [12.0, 10.0, 8.0, 6.0, 4.0])))

    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)

    in_box = np.zeros((H, W), dtype=bool)
    in_box[int(r0) + 4:int(r1) - 4, int(c0) + 4:W - 4] = True
    # Exclude top-left legend
    in_box[:int(r0 + 0.22 * (r1 - r0)), :int(c0 + 0.35 * (W - c0))] = False

    rline = in_box & (R > 165.0) & (G < 85.0) & (B < 85.0)
    bline = in_box & (B > 125.0) & (R < 100.0) & (B - R > 65.0)
    # Any shaded ribbon pixel (pink, pale blue, or overlapping purplish blend) is tinted below pure white
    nonwhite = in_box & ((R < 248.0) | (G < 248.0) | (B < 248.0))
    rrib_pure = nonwhite & (R - B > 8.0) & (R > 165.0) & ~bline
    brib_pure = nonwhite & (B - R > 4.0) & (B > 165.0) & ~rline

    line_mask = rline if series == "control" else bline
    other_line = bline if series == "control" else rline
    pure_rib = rrib_pure if series == "control" else brib_pure

    out = []
    for xv in at:
        col_f = xa.to_pixel(float(xv))
        col = int(round(col_f))
        # Find curve line center (slightly wider search at endpoint 700 nm where curve ends ~2 px left)
        rows = np.array([], dtype=int)
        for hw in (1, 3, 5):
            c_lo = max(int(c0) + 4, col - hw)
            c_hi = min(W - 4, col + hw + 1)
            rows = np.where(line_mask[:, c_lo:c_hi].any(axis=1))[0]
            if len(rows) > 0:
                break
        if len(rows) == 0:
            # At exact line crossing (e.g. ~400 nm or ~680 nm), blue may overpaint red: check adjacent cols +-8 px
            c_lo = max(int(c0) + 4, col - 8)
            c_hi = min(W - 4, col + 9)
            rows = np.where(line_mask[:, c_lo:c_hi].any(axis=1))[0]
        if len(rows) == 0:
            out.append({"x": float(xv), "y": None, "note": f"no {series} curve at col {col}"})
            continue
        row_c = float(rows.mean() + 0.5)
        half_px = max(0.5 * (rows.max() - rows.min() + 1), 2.0)

        # Find symmetric SE ribbon around row_c (at endpoints 350/700 nm step inward up to 6 px so ribbon is hit)
        col_rib = max(int(c0) + 4, min(W - 5, col))
        if not nonwhite[:, col_rib].any():
            for dcol in (-3, 3, -6, 6, -10, 10):
                cc = max(int(c0) + 4, min(W - 5, col + dcol))
                if nonwhite[:, cc].any():
                    col_rib = cc
                    break
        elif xv >= 695.0:
            col_rib = max(int(c0) + 4, col - 6)
        elif xv <= 355.0:
            col_rib = min(W - 5, col + 6)
        pure_col = pure_rib[:, col_rib] | line_mask[:, col_rib]
        nw_col = nonwhite[:, col_rib]

        def span_around(mask_1d: np.ndarray, rc: float, max_gap: int = 8) -> tuple[float, float] | None:
            idx = np.where(mask_1d)[0]
            if len(idx) == 0:
                return None
            for s, e, _ in find_runs(idx, gap=max_gap):
                if s - 10 <= rc <= e + 10:
                    return (float(s + 0.5), float(e + 0.5))
            return None

        sp_pure = span_around(pure_col, row_c)
        sp_nw = span_around(nw_col, row_c)
        # Check if the other series is nearby (overlapping ribbons)
        other_rows = np.where(other_line[:, max(int(c0) + 4, col_rib - 5):min(W - 4, col_rib + 6)].any(axis=1))[0]
        other_rc = float(other_rows.mean() + 0.5) if len(other_rows) else None

        # Determine if red and blue ribbons overlap at this wavelength:
        # Compare pure_rib span vs full nonwhite span; if nonwhite extends much further toward other_rc,
        # the two ribbons overlap and the outer edge of nonwhite AWAY from other_rc is our unoccluded SE edge.
        se_px = None
        sym_note = ""
        if sp_pure is not None and sp_nw is not None and other_rc is not None:
            pure_h = sp_pure[1] - sp_pure[0]
            nw_h = sp_nw[1] - sp_nw[0]
            overlap = (nw_h > pure_h + 20.0) or (abs(other_rc - row_c) < 280.0)
            if overlap:
                if row_c <= other_rc:
                    # Our curve is ABOVE (smaller row) -> upper edge of nonwhite union is unoccluded
                    se_px = max(row_c - sp_nw[0], 10.0)
                else:
                    # Our curve is BELOW (larger row) -> lower edge of nonwhite union is unoccluded
                    se_px = max(sp_nw[1] - row_c, 10.0)
                r_top, r_bot = row_c - se_px, row_c + se_px
                sym_note = " (symmetric mirror across overlap)"
        if se_px is None:
            if sp_pure is not None:
                r_top, r_bot = sp_pure
            elif sp_nw is not None:
                r_top, r_bot = sp_nw
                sym_note = " (blended overlap)"
            else:
                r_top, r_bot = row_c - half_px, row_c + half_px
                sym_note = " (no ribbon)"

        y_val, lo, hi, _, _ = y_interval(ya, row_c, half_px)
        bar_hi = float(ya.to_value(r_top))
        bar_lo = float(ya.to_value(r_bot))
        h_rib = 4.0 + max(ya.residual_px, 0.5)
        bar_hi_outer = float(ya.to_value(r_top - h_rib))
        bar_lo_outer = float(ya.to_value(r_bot + h_rib))
        out.append({
            "x": float(xv),
            "y": float(y_val), "lo": float(lo), "hi": float(hi),
            "bar_lo": float(min(bar_lo, bar_hi)), "bar_hi": float(max(bar_lo, bar_hi)),
            "bar_lo_outer": float(min(bar_lo_outer, bar_hi_outer)),
            "bar_hi_outer": float(max(bar_lo_outer, bar_hi_outer)),
            "px": [float(col_f), row_c],
            "px_lo": float(r_top), "px_hi": float(r_bot),
            "note": f"{series} mean +- SE{sym_note} (rx={xa.residual_px:.2f} px, ry={ya.residual_px:.2f} px)",
        })
    return out


def _calibrate_8no(image: np.ndarray, y_labels: list[float]) -> tuple[D.Axis, float, float, list[float]]:
    gray = image.astype(float).mean(axis=2)
    dark = gray < 100.0
    c_frac = dark.mean(axis=0)
    c_frac[:100] = 0.0
    c0 = run_centers(np.where(c_frac > 0.50)[0], gap=5)[0]
    r1 = run_centers(np.where(dark.mean(axis=1) > 0.50)[0], gap=5)[-1]
    yt = run_centers(np.where(dark[:int(r1), int(round(c0)) - 20])[0], gap=4)
    xt = run_centers(np.where(dark[int(round(r1)) + 20, int(round(c0)):])[0] + int(round(c0)), gap=4)
    if len(yt) != len(y_labels) or len(xt) != 2:
        raise RuntimeError(f"Fig. 8N/O expected {len(y_labels)} y-ticks and 2 x-ticks, got yt={yt}, xt={xt}")
    ya = D.Axis("linear", list(zip(yt, y_labels)))
    return ya, c0, r1, xt


def read_8no_mean_sd(image: np.ndarray, y_labels: list[float], what: str) -> list[dict]:
    """Read red mean dot and vertical red SD line for control (x=1) and cyto (x=2) in Fig. 8N or 8O."""
    ya, c0, r1, xt = _calibrate_8no(image, y_labels)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    red = (R - G > 80.0) & (R - B > 80.0)

    out = []
    for idx, (xc, lbl) in enumerate(zip(xt, ["control", "cyto"]), start=1):
        sub_red = red[:int(r1) - 10, int(round(xc)) - 80:int(round(xc)) + 80]
        rows_red = np.where(sub_red.any(axis=1))[0]
        widths = sub_red.sum(axis=1)
        dot_rows = np.where(widths > widths[rows_red].min() + 15)[0]
        row_mean = float(dot_rows.mean() + 0.5)
        r_top = float(rows_red.min() + 0.5)
        r_bot = float(rows_red.max() + 0.5)

        y_mean, lo, hi, _, _ = y_interval(ya, row_mean, 2.0)
        bar_hi = float(ya.to_value(r_top))
        bar_lo = float(ya.to_value(r_bot))
        h_px = 2.0 + max(ya.residual_px, 0.5)
        bar_hi_outer = float(ya.to_value(r_top - h_px))
        bar_lo_outer = float(ya.to_value(r_bot + h_px))
        sd_val = 0.5 * abs(bar_hi - bar_lo)
        out.append({
            "x": float(idx),
            "y": float(y_mean), "lo": float(lo), "hi": float(hi),
            "bar_lo": float(min(bar_lo, bar_hi)), "bar_hi": float(max(bar_lo, bar_hi)),
            "bar_lo_outer": float(min(bar_lo_outer, bar_hi_outer)),
            "bar_hi_outer": float(max(bar_lo_outer, bar_hi_outer)),
            "px": [float(xc), row_mean],
            "px_lo": float(r_top), "px_hi": float(r_bot),
            "note": f"{lbl} {what} mean={y_mean:.5g}, SD={sd_val:.5g} (ry={ya.residual_px:.2f} px)",
        })
    return out


def read_8no_dots(image: np.ndarray, y_labels: list[float], group_index: int,
                  what: str) -> list[dict]:
    """Detect all individual stacked circular black dots for control (1) or cyto (2) in Fig. 8N or 8O
    via smoothed Euclidean distance transform with greedy non-maximum suppression."""
    ya, c0, r1, xt = _calibrate_8no(image, y_labels)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    # Unified seamless mask of black dots + overlaid red SD stem/dot (eliminates anti-aliased seam)
    ink = (G < 100.0) & (B < 100.0)
    blk = (R < 40.0) & (G < 40.0) & (B < 40.0)

    # Measure isolated dot diameter directly from topmost circular black blob inside plot
    lab0, n0 = ndimage.label(blk[:int(r1) - 15, int(c0) + 40:])
    top_r, diam = 99999, 135.0
    for k in range(1, n0 + 1):
        rr, cc = np.where(lab0 == k)
        if len(rr) > 3000:
            h, w = rr.max() - rr.min() + 1, cc.max() - cc.min() + 1
            if 0.85 < h / w < 1.15 and rr.min() < top_r:
                top_r = int(rr.min())
                diam = 0.5 * (h + w)
    rad = diam / 2.0

    # Find the widest empty white column gap between xt[0] and xt[1] to separate groups cleanly
    c_mid_lo, c_mid_hi = int(round(xt[0] + 150)), int(round(xt[1] - 150))
    col_has_ink = ink[:int(r1) - 15, c_mid_lo:c_mid_hi].any(axis=0)
    empty_runs = find_runs(np.where(~col_has_ink)[0], gap=1)
    if empty_runs:
        best_gap = max(empty_runs, key=lambda r: r[1] - r[0])
        c_split = int(round(best_gap[2])) + c_mid_lo
    else:
        c_split = int(round(0.5 * (xt[0] + xt[1])))

    if group_index == 1:
        cl, ch = int(c0) + 20, c_split
    else:
        cl, ch = c_split, image.shape[1] - 10
    sub = ink[:int(r1) - 15, cl:ch]
    dist = ndimage.gaussian_filter(ndimage.distance_transform_edt(sub), sigma=4.0)
    fp = int(round(rad * 1.10))
    if fp % 2 == 0:
        fp += 1
    mx = ndimage.maximum_filter(dist, size=fp)
    pks = (dist == mx) & (dist > 0.72 * rad)
    plab, pn = ndimage.label(pks)

    cands = []
    for k in range(1, pn + 1):
        rr, cc = np.where(plab == k)
        rc = float(rr.mean() + 0.5)
        cc_m = float(cc.mean() + 0.5) + cl
        cands.append((float(dist[rr, cc].max()), rc, cc_m))

    # Greedy NMS by peak distance, rejecting saddle points closer than 1.05 * rad to a stronger dot center
    cands.sort(key=lambda t: -t[0])
    kept = []
    min_sep = 1.05 * rad
    for dval, rc, cc_m in cands:
        if all((rc - r2) ** 2 + (cc_m - c2) ** 2 >= min_sep ** 2 for _, r2, c2 in kept):
            kept.append((dval, rc, cc_m))

    # Sort ascending by y value (descending pixel row)
    kept.sort(key=lambda t: (-t[1], t[2]))
    lbl = "control" if group_index == 1 else "cyto"
    out = []
    for idx, (_, rc, cc_m) in enumerate(kept, start=1):
        y_val, lo, hi, px_lo, px_hi = y_interval(ya, rc, 3.0)
        out.append({
            "x": float(idx), "y": float(y_val), "lo": float(lo), "hi": float(hi),
            "px": [float(cc_m), float(rc)], "px_lo": px_lo, "px_hi": px_hi,
            "note": f"{lbl} {what} dot {idx}/{len(kept)} (ry={ya.residual_px:.2f} px)",
        })
    return out
