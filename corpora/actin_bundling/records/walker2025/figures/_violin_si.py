"""Shared calibrated reader for Fig. S7F and S7G log-scale multimer length violin plots (900 dpi).

- Calibrates each stacked panel's log10 y-axis directly from the 3 long major outward y-ticks
  (10^2, 10^1, 10^0) which extend beyond the shorter minor ticks (detected by probing just past
  the minor-tick length without rejecting ticks that touch their numeric labels).
- Fig. S7G (mode="white_dot"): locates seaborn inner="box" white median circle on each violin center
  column (or collapsed flat baseline at y = 1.0 when median = 1).
- Fig. S7F (mode="quartile_dash"): locates seaborn inner="quartile" horizontal black dashed median
  line inside the colored violin interior using relative luminance drop below the violin's fill gray
  (recording Q1/Q3 in note; returns 1.0 when collapsed flat at y = 10^0).
"""

from __future__ import annotations

import numpy as np
import _common as C
from _bars_stacked import find_stacked_boxes


def _find_log_major_ticks(img: np.ndarray, box: dict, n_decades: int = 3) -> list[float]:
    """Return pixel row centers of the `n_decades` long major log ticks (top-to-bottom: 10^2, 10^1, 10^0)."""
    dark = img.sum(axis=2, dtype=np.uint16) < 345
    c_left_lo = box["left"][0]
    r_lo = max(0, box["top"]["lo"] - 6)
    r_hi = min(dark.shape[0], box["bot"]["hi"] + 7)
    strip = dark[r_lo:r_hi, max(0, c_left_lo - 45) : c_left_lo][:, ::-1]
    # Columns 26..38 lie strictly past minor ticks (<=24 px) and right inside the long major ticks!
    deep = strip[:, 26:38].mean(axis=1) > 0.5
    cls = [c[0] + r_lo for c in C.clusters_1d(np.flatnonzero(deep), gap=3) if (c[2] - c[1] + 1) <= 22]
    return C.match_regular_ticks(cls, n_decades, tol=8.0)


def read_violin_panel(
    img: np.ndarray,
    n_rows: int,
    row_idx: int,
    mode: str,
    x_vals: list[float] = (10.0, 1.0, 0.1, 0.01, 0.001),
    label: str = "",
) -> list[dict]:
    boxes = find_stacked_boxes(img, n_rows, max_tick_len=35)
    b = boxes[row_idx]
    c0, r0, c1, r1 = b["box"]

    yt = _find_log_major_ticks(img, b, n_decades=3)
    ya = C.make_axis("log10", list(zip(yt, [100.0, 10.0, 1.0])), max_res_px=4.5)

    xt = C.match_regular_ticks(boxes[-1]["x_ticks"], len(x_vals), tol=6.0)
    col_sp = (xt[-1] - xt[0]) / max(len(x_vals) - 1, 1)
    half_w = int(round(0.42 * col_sp))
    r_one = float(ya.to_pixel(1.0))

    rgb = img.astype(np.float64)
    gray = rgb.mean(axis=2)
    sat = rgb.max(axis=2) - rgb.min(axis=2)

    out = []
    for vi, xv in enumerate(x_vals):
        cx = xt[vi]
        ci = int(round(cx))
        vc0 = max(c0 + 2, ci - half_w)
        vc1 = min(c1 - 2, ci + half_w)

        if mode == "white_dot":
            # Fig. S7G: crisp bright white circle dead-center on violin column [ci - 5 : ci + 6],
            # tightly enclosed on left, right, up, and down by black boxplot/outline ink!
            # Restrict vertical search to y >= 0.85 (row <= ya.to_pixel(0.85)) so we never hit bottom spine.
            r_search_hi = min(r1 - 2, int(round(r_one + 14)))
            strip_g = gray[r0 + 2 : r_search_hi, ci - 4 : ci + 5]
            strip_s = sat[r0 + 2 : r_search_hi, ci - 4 : ci + 5]
            is_white = (strip_g > 200.0) & (strip_s < 45.0)
            # Flanked tightly on BOTH left and right within 11 px by dark boxplot wall
            is_dark_l = (gray[r0 + 2 : r_search_hi, ci - 11 : ci - 1] < 110.0).any(axis=1)
            is_dark_r = (gray[r0 + 2 : r_search_hi, ci + 2 : ci + 12] < 110.0).any(axis=1)
            white_rows = np.flatnonzero(is_white.any(axis=1) & is_dark_l & is_dark_r)
            scored = []
            for ctr, s_r, e_r in C.clusters_1d(white_rows, gap=1):
                h_dot = e_r - s_r + 1
                if 1 <= h_dot <= 14:
                    r_abs = ctr + (r0 + 2)
                    rr_i = int(round(r_abs))
                    u_blk = (gray[max(r0, rr_i - 10) : max(r0 + 1, rr_i - 2), ci - 3 : ci + 4] < 110.0).any()
                    d_blk = (gray[min(r1 - 1, rr_i + 3) : min(r1, rr_i + 11), ci - 3 : ci + 4] < 110.0).any()
                    # Check width of white dot along center row is small (<= 12 px, not wide background white!)
                    w_white = (gray[rr_i, ci - 12 : ci + 13] > 200.0).sum()
                    if (u_blk or d_blk) and w_white <= 12:
                        peak_g = float(gray[max(r0, rr_i - 2) : min(r1, rr_i + 3), ci - 3 : ci + 4].max())
                        score = int(u_blk) + int(d_blk)
                        scored.append((score, peak_g, -abs(h_dot - 5), r_abs))
            if scored:
                scored.sort(reverse=True)
                row_med = scored[0][3]
                note = f"{label} row={row_idx} col={vi+1}: white-dot median, ya_res={ya.residual_px:.2f}px"
            else:
                row_med = r_one
                note = f"{label} row={row_idx} col={vi+1}: collapsed at y=1.0, ya_res={ya.residual_px:.2f}px"
            out.append(
                C.reading_point(
                    x_val=float(xv),
                    col=cx,
                    row=row_med,
                    y_axis=ya,
                    half_px=3.0,
                    note=note,
                    clamp_min=1.0,
                )
            )

        elif mode == "quartile_dash":
            # Fig. S7F: seaborn inner="quartile" horizontal dashed black lines inside colored violin body.
            sub_g = gray[r0 + 3 : r1 - 2, vc0 : vc1 + 1]
            sub_s = sat[r0 + 3 : r1 - 2, vc0 : vc1 + 1]
            colored = (sub_s > 30.0) & (sub_g < 220.0)
            if colored.sum() < 400:
                out.append(
                    C.reading_point(
                        x_val=float(xv), col=cx, row=r_one, y_axis=ya, half_px=3.5,
                        note=f"{label} row={row_idx} col={vi+1}: collapsed at y=1.0, ya_res={ya.residual_px:.2f}px",
                        clamp_min=1.0,
                    )
                )
                continue

            # Measure fill luminance of this violin color; dashed lines are significantly darker than fill!
            fill_gray = float(np.percentile(sub_g[colored], 65))
            dark_thresh = min(110.0, fill_gray - 18.0)
            dark_rel = sub_g < dark_thresh

            dash_rows = []
            for ri_loc in range(sub_g.shape[0]):
                c_cols = np.flatnonzero(colored[ri_loc])
                if c_cols.size < 18:
                    continue
                span = int(c_cols[-1]) - int(c_cols[0])
                if span < 22:
                    continue
                # Scan strictly inside the middle ~72% of the colored violin body (immune to dark left/right borders!)
                pad = max(6, int(round(0.14 * span)))
                il = int(c_cols[0]) + pad
                ir = int(c_cols[-1]) - pad
                if ir - il < 10:
                    continue
                row_d = dark_rel[ri_loc, il : ir + 1]
                frac_d = float(row_d.mean())
                if 0.14 <= frac_d <= 0.90:
                    runs = C.clusters_1d(np.flatnonzero(row_d), gap=2)
                    if not runs:
                        continue
                    avg_run = float(np.mean([t[2] - t[1] + 1 for t in runs]))
                    max_run = float(max(t[2] - t[1] + 1 for t in runs))
                    single_interior_long = (
                        len(runs) == 1
                        and 18.0 <= max_run <= 38.0
                        and runs[0][1] >= 4
                        and runs[0][2] <= (ir - il) - 4
                    )
                    if len(runs) >= 2 or single_interior_long or (ir - il <= 28 and len(runs) >= 1 and frac_d >= 0.30):
                        dash_rows.append((ri_loc + (r0 + 3), frac_d, avg_run, max_run, len(runs)))

            # Cluster consecutive dash rows into distinct horizontal quartile lines
            lines = []
            if dash_rows:
                groups = []
                for item in dash_rows:
                    if not groups or item[0] > groups[-1][-1][0] + 6:
                        groups.append([item])
                    else:
                        groups[-1].append(item)
                for g_i in groups:
                    if g_i[-1][0] - g_i[0][0] <= 12:
                        r_ctr = 0.5 * (g_i[0][0] + g_i[-1][0]) + 0.5
                        avg_dash = max(x[2] for x in g_i)
                        max_dash = max(x[3] for x in g_i)
                        lines.append((r_ctr, avg_dash, max_dash, float(ya.to_value(r_ctr))))

            # Keep lines cleanly above y = 1.05 (r < r_one - 4)
            above = [ln for ln in lines if ln[0] < r_one - 4.0]
            above.sort(key=lambda t: t[0])  # top-to-bottom (descending y)
            if len(above) >= 3:
                if len(above) == 3:
                    q3_r, med_r, q1_r = above[0][0], above[1][0], above[2][0]
                else:
                    interior = above[1:-1]
                    med_ln = max(interior, key=lambda t: (t[2], t[1]))
                    med_r = med_ln[0]
                    q3_r, q1_r = above[0][0], above[-1][0]
                row_med = med_r
                q_str = f"Q1={float(ya.to_value(q1_r)):.2f}, Q3={float(ya.to_value(q3_r)):.2f}"
            elif len(above) == 2:
                # If lower line is >= 2.5, it is Q1 and upper line is Median (with Q3 in the thin upper spire)
                if above[1][3] >= 2.5 or above[0][2] > above[1][2] * 1.05:
                    row_med = above[0][0]
                    q_str = f"Q1={above[1][3]:.2f}, Med={above[0][3]:.2f}"
                else:
                    row_med = above[1][0]
                    q_str = f"Q1~1.0, Q3={above[0][3]:.2f}"
            elif len(above) == 1:
                if above[0][2] >= 12.0:
                    row_med = above[0][0]
                    q_str = f"Q1~1.0, Med={above[0][3]:.2f}"
                else:
                    row_med = r_one
                    q_str = f"Q1=1.0, Med=1.0, Q3={above[0][3]:.2f}"
            else:
                row_med = r_one
                q_str = "collapsed at y=1.0"

            out.append(
                C.reading_point(
                    x_val=float(xv),
                    col=cx,
                    row=row_med,
                    y_axis=ya,
                    half_px=3.5,
                    note=f"{label} row={row_idx} col={vi+1}: {q_str}, ya_res={ya.residual_px:.2f}px",
                    clamp_min=1.0,
                )
            )
        else:
            raise ValueError(f"Unknown mode: {mode}")
    return out
