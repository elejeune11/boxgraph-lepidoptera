"""Calibrated reader for Fig. 2E, 2F, 2G biomimetic bead assay bar charts."""

from __future__ import annotations

import numpy as np

import _calib as C


def _read_two_bars(
    image: np.ndarray,
    tick_vals: list[float],
    tick_min_step: float,
    labels: list[str],
    panel_name: str,
    bar2_has_err: bool = True,
) -> list[dict]:
    dark, gray, h_lines, v_lines = C.find_l_spines(image, thresh=140, row_frac=0.30, col_frac=0.30)
    x_spine = max(h_lines, key=lambda t: t[0])
    y_spine = min(v_lines, key=lambda t: t[0])

    yt_raw = C.find_ticks(dark, y_spine, (0, int(x_spine[1]) - 50), "y", side=+1, band=(4, 24), frac=0.35, max_width=50)
    yt = C.match_ticks_regular(yt_raw, len(tick_vals), tol=6.0, min_step=tick_min_step)
    y_axis = C.verify_axis("linear", list(zip(yt, tick_vals)), max_residual_px=2.5)

    H, W = gray.shape
    c0 = int(y_spine[2]) + 15
    r_base = int(x_spine[1]) - 15
    band = gray[r_base - 80 : r_base, c0:W] < 235
    col_fill = band.mean(axis=0) > 0.7
    runs = [r for r in C.runs_1d(col_fill, gap=3) if (r[1] - r[0] + 1) >= 150]
    if len(runs) != 2:
        raise RuntimeError(f"{panel_name}: expected 2 bars, found {runs}")

    out = []
    for idx, (a, b) in enumerate(runs, start=1):
        ca, cb = a + c0, b + c0
        c_mid = 0.5 * (ca + cb) + 0.5
        margin = int(0.12 * (cb - ca))
        left_strip = gray[:r_base, ca + margin : int(c_mid) - 35]
        right_strip = gray[:r_base, int(c_mid) + 35 : cb - margin]
        # Scan UPWARD from r_base - 1 so text labels above the bar (e.g. 'Comet Tail' in Fig 2G) are never hit!
        is_bar_row = (left_strip.mean(axis=1) < 230) & (right_strip.mean(axis=1) < 230)
        top_r = r_base - 1
        while top_r > 0 and is_bar_row[top_r - 1]:
            top_r -= 1
        top_row = float(top_r)

        has_err = (idx == 1) or bar2_has_err
        u_row, l_row = None, None
        if has_err:
            err_mask_up = gray < 120
            err_mask_up[int(top_row) + 1 :, :] = False
            u_row, _ = C.find_vertical_error_bar(
                err_mask_up,
                c_mid,
                top_row,
                row_min=5,
                row_max=int(top_row),
                stem_half_w=12,
                max_gap=8,
                cap_half_w=35,
                min_cap_width=18,
                start_up_offset=2,
                start_dn_offset=0,
            )
            bar_gray = float(np.median(left_strip[min(r_base - 10, int(top_row) + 30) : min(r_base, int(top_row) + 120)]))
            err_mask_dn = gray < (0.55 * bar_gray)
            err_mask_dn[: int(top_row), :] = False
            _, l_row = C.find_vertical_error_bar(
                err_mask_dn,
                c_mid,
                top_row,
                row_min=int(top_row),
                row_max=r_base,
                stem_half_w=12,
                max_gap=8,
                cap_half_w=35,
                min_cap_width=18,
                start_up_offset=0,
                start_dn_offset=2,
            )
            if u_row is not None and (top_row - u_row) < 6:
                u_row = None
            if l_row is not None and (l_row - top_row) < 6:
                l_row = None

        note = f"{panel_name} bar {idx} ({labels[idx - 1]}); cal_res_y={y_axis.residual_px:.2f}px"
        out.append(
            C.make_point_reading(
                float(idx),
                c_mid,
                top_row,
                y_axis,
                upper_err_row=u_row,
                lower_err_row=l_row,
                mirror_lower=True,
                err_px=2.0 + y_axis.residual_px,
                note=note,
            )
        )
    return out


def read_panel(image: np.ndarray, panel: str) -> list[dict]:
    if panel == "2e":
        return _read_two_bars(
            image,
            [3.0, 2.5, 2.0, 1.5, 1.0, 0.5],
            220.0,
            ["-alpha-Actinin", "+alpha-Actinin"],
            "Fig. 2E",
            bar2_has_err=True,
        )
    if panel == "2f":
        return _read_two_bars(
            image,
            [1.0, 0.8, 0.6, 0.4, 0.2],
            280.0,
            ["-alpha-Actinin", "+alpha-Actinin"],
            "Fig. 2F",
            bar2_has_err=True,
        )
    if panel == "2g":
        return _read_two_bars(
            image,
            [1.0, 0.8, 0.6, 0.4, 0.2],
            260.0,
            ["Comet Tail", "Protrusions (normalized to 1.0, no error bar)"],
            "Fig. 2G",
            bar2_has_err=False,
        )
    raise ValueError(f"Unknown panel {panel!r}")
