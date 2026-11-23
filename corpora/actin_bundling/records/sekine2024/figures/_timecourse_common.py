"""Shared calibration and line+ribbon reading helpers for time-course plots (Fig. 1h-k, Supp. Fig. 1g)."""

from __future__ import annotations

import numpy as np


def _split_runs(idx: np.ndarray, gap: int = 5) -> list[np.ndarray]:
    if len(idx) == 0:
        return []
    cuts = np.where(np.diff(idx) > gap)[0] + 1
    return [g for g in np.split(idx, cuts) if len(g) > 0]


def _regular_ticks(cands: list[float], tol: float = 5.0, min_gap: float = 100.0) -> list[float]:
    c = sorted(cands)
    if len(c) <= 2:
        return c
    best: list[float] = []
    for i in range(len(c)):
        for j in range(i + 1, len(c)):
            step = c[j] - c[i]
            if step < min_gap:
                continue
            members = [c[i], c[j]]
            while True:
                est = (members[-1] - members[0]) / (len(members) - 1)
                target = members[-1] + est
                within = max(tol, 0.02 * est)
                near = [x for x in c if abs(x - target) <= within]
                if not near:
                    break
                members.append(min(near, key=lambda x: abs(x - target)))
            if len(members) > len(best) or (
                len(members) == len(best) and len(best) > 1 and (members[1] - members[0] > best[1] - best[0])
            ):
                best = members
    return best if len(best) >= 2 else c


def calibrate_l_axes(
    image: np.ndarray,
    x_labels: list[float],
    y_labels_top_to_bottom: list[float],
    tick_band_x: tuple[int, int] = (2, 11),
    tick_band_y: tuple[int, int] = (2, 11),
    min_gap_x: float = 120.0,
    min_gap_y: float = 120.0,
    dark_thresh: float = 150.0,
    tick_thresh: float = 25.0,
):
    """Calibrate an L-shaped (bottom x, left y) axis pair from outward ticks."""
    img = image.astype(float)
    H, W, _ = img.shape
    gray = img.mean(axis=2)
    dark = gray < dark_thresh

    hr_groups = _split_runs(np.where(dark.sum(axis=1) > W * 0.40)[0], gap=10)
    vc_groups = _split_runs(np.where(dark.sum(axis=0) > H * 0.40)[0], gap=10)
    if not hr_groups or not vc_groups:
        raise RuntimeError(f"Could not find L-axes (hr={len(hr_groups)}, vc={len(vc_groups)})")

    hr = hr_groups[-1]  # bottom x-axis
    vc = vc_groups[0]   # left y-axis

    r0_xt = int(hr.max()) + tick_band_x[0]
    r1_xt = min(H, int(hr.max()) + tick_band_x[1])
    xp = (255.0 - gray[r0_xt:r1_xt, :]).mean(axis=0)
    xp[: int(vc.max()) + 2] = 0.0
    xt_all = [float(((g + 0.5) * xp[g]).sum() / xp[g].sum()) for g in _split_runs(np.where(xp > tick_thresh)[0], gap=4)]
    xt = _regular_ticks(xt_all, tol=5.0, min_gap=min_gap_x)
    if len(xt) != len(x_labels):
        raise RuntimeError(f"Expected {len(x_labels)} x-ticks, found {len(xt)} from {xt_all}")

    c1_yt = max(0, int(vc.min()) - tick_band_y[0])
    c0_yt = max(0, int(vc.min()) - tick_band_y[1])
    yp = (255.0 - gray[:, c0_yt:c1_yt]).mean(axis=1)
    yp[int(hr.min()) + 12 :] = 0.0
    yt_all = [float(((g + 0.5) * yp[g]).sum() / yp[g].sum()) for g in _split_runs(np.where(yp > tick_thresh)[0], gap=4)]
    yt = _regular_ticks(yt_all, tol=5.0, min_gap=min_gap_y)
    if len(yt) != len(y_labels_top_to_bottom):
        raise RuntimeError(f"Expected {len(y_labels_top_to_bottom)} y-ticks, found {len(yt)} from {yt_all}")

    bx, ax = np.polyfit(xt, x_labels, 1)
    by, ay = np.polyfit(yt, y_labels_top_to_bottom, 1)
    res_x = float(np.max(np.abs(np.array(xt) - (np.array(x_labels) - ax) / bx)))
    res_y = float(np.max(np.abs(np.array(yt) - (np.array(y_labels_top_to_bottom) - ay) / by)))
    if res_x > 2.5 or res_y > 2.5:
        raise RuntimeError(f"Calibration residual too high: rx={res_x:.2f}px, ry={res_y:.2f}px")

    return {
        "box": (int(vc.max()) + 4, 20, W - 20, int(hr.min()) - 4),
        "ax": float(ax),
        "bx": float(bx),
        "ay": float(ay),
        "by": float(by),
        "res_x": res_x,
        "res_y": res_y,
    }


def read_timecourse_ribbon(
    image: np.ndarray,
    at: list[float],
    x_labels: list[float],
    y_labels_top_to_bottom: list[float],
    color_kind: str = "green",
    tick_band_x: tuple[int, int] = (2, 11),
    tick_band_y: tuple[int, int] = (2, 11),
    min_gap_x: float = 120.0,
    min_gap_y: float = 120.0,
):
    """Read mean line + shaded SD/range ribbon at each requested time point."""
    cal = calibrate_l_axes(
        image, x_labels, y_labels_top_to_bottom,
        tick_band_x=tick_band_x, tick_band_y=tick_band_y,
        min_gap_x=min_gap_x, min_gap_y=min_gap_y,
    )
    c0, r0, c1, r1 = cal["box"]
    ax, bx, ay, by = cal["ax"], cal["bx"], cal["ay"], cal["by"]

    sub = image[r0:r1, c0:c1].astype(float)
    R, G, B = sub[:, :, 0], sub[:, :, 1], sub[:, :, 2]

    if color_kind == "green":
        line_mask = (G - R > 45.0) & (G - B > 45.0) & (G > 70.0)
        ribbon_mask = ((G - R >= 4.0) & (G - B >= 4.0) & (G > 195.0)) | line_mask
    elif color_kind == "blue":
        line_mask = (B - R > 80.0) & (B - G > 80.0) & (B > 120.0)
        ribbon_mask = ((B - R >= 5.0) & (B - G >= 5.0) & (B > 210.0)) | line_mask
    else:
        raise ValueError(f"Unknown color_kind: {color_kind}")

    W_s = sub.shape[1]
    out = []
    xmin, xmax = min(at), max(at)
    for xv in at:
        xv = float(xv)
        col_f = (xv - ax) / bx
        col_s = int(round(col_f)) - c0
        # Slight inward nudge only at endpoints so line/ribbon edge is never missed
        if abs(xv - xmin) < 1e-6:
            c_line_rng = range(max(0, col_s), min(W_s, col_s + 5))
            c_rib = min(W_s - 1, max(0, col_s + 3))
        elif abs(xv - xmax) < 1e-6:
            c_line_rng = range(max(0, col_s - 4), min(W_s, col_s + 1))
            c_rib = min(W_s - 1, max(0, col_s - 3))
        else:
            c_line_rng = range(max(0, col_s - 3), min(W_s, col_s + 4))
            c_rib = min(W_s - 1, max(0, col_s))

        # Mean line center (median of column midpoints over narrow window around col_s)
        mids = []
        spans = []
        for cc in c_line_rng:
            rr = np.nonzero(line_mask[:, cc])[0]
            if len(rr):
                runs = _split_runs(rr, gap=40)  # bridges hollow circle center (~20px)
                best = max(runs, key=lambda g: int(g[-1] - g[0] + 1))
                mids.append(0.5 * (float(best[0]) + float(best[-1])) + 0.5)
                spans.append(int(best[-1] - best[0] + 1))
        if not mids:
            out.append({"x": xv, "y": None, "note": "curve not found"})
            continue

        rc_abs = float(np.median(mids)) + r0
        y_val = ay + by * rc_abs
        err_px = 2.0 + cal["res_y"]
        err_y = abs(by) * err_px

        # Ribbon upper/lower edges (use 3 nearby columns to skip black vertical dotted line at x=0)
        rib_cols = [c for c in (c_rib - 3, c_rib, c_rib + 3) if 0 <= c < W_s]
        if abs(xv - xmin) < 1e-6:
            rib_cols = [c for c in (c_rib, c_rib + 2, c_rib + 4) if 0 <= c < W_s]
        elif abs(xv - xmax) < 1e-6:
            rib_cols = [c for c in (c_rib - 4, c_rib - 2, c_rib) if 0 <= c < W_s]

        top_rows, bot_rows = [], []
        for cc in rib_cols:
            rr = np.nonzero(ribbon_mask[:, cc])[0]
            if len(rr):
                runs = _split_runs(rr, gap=45)
                # pick the run containing rc_abs - r0
                containing = [g for g in runs if g[0] - 5 <= (rc_abs - r0) <= g[-1] + 5]
                best = containing[0] if containing else max(runs, key=len)
                top_rows.append(float(best[0]) + r0)
                bot_rows.append(float(best[-1]) + 1.0 + r0)

        row_dict = {
            "x": xv,
            "y": y_val,
            "lo": y_val - err_y,
            "hi": y_val + err_y,
            "px": [col_f, rc_abs],
            "px_lo": rc_abs - err_px,
            "px_hi": rc_abs + err_px,
            "note": f"{color_kind} line+ribbon (rx={cal['res_x']:.2f}px, ry={cal['res_y']:.2f}px)",
        }
        if top_rows and bot_rows:
            r_hi_val = float(np.median(top_rows))  # smaller row = higher y value
            r_lo_val = float(np.median(bot_rows))  # larger row = lower y value
            bar_hi = ay + by * r_hi_val
            bar_lo = ay + by * r_lo_val
            bar_pad = abs(by) * (3.0 + cal["res_y"])
            row_dict["bar_lo"] = min(bar_lo, y_val)
            row_dict["bar_hi"] = max(bar_hi, y_val)
            row_dict["bar_lo_outer"] = row_dict["bar_lo"] - bar_pad
            row_dict["bar_hi_outer"] = row_dict["bar_hi"] + bar_pad
            row_dict["px_lo"] = r_hi_val
            row_dict["px_hi"] = r_lo_val
        out.append(row_dict)
    return out
