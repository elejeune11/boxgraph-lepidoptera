"""Shared calibrated reader for Fig. 6D experimental stacked bar chart with SD error bars (1200 dpi)."""

from __future__ import annotations
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _plot_common as pc

RATIOS = [20.0, 10.0, 7.5]


def _find_error_bar(eb_mask: np.ndarray, r0: int, anchor_row: int, x_lo: int, x_hi: int):
    """Find the vertical error bar stem crossing within +-6 px of anchor_row in columns [x_lo, x_hi)."""
    best = (0, None, None, None)
    ar = anchor_row - r0
    H = eb_mask.shape[0]
    for cc in range(x_lo, min(x_hi, eb_mask.shape[1])):
        win = np.flatnonzero(eb_mask[max(0, ar - 6):min(H, ar + 7), cc])
        if not win.size:
            continue
        r_s = max(0, ar - 6) + int(win[0])
        ru, rd = r_s, r_s
        while ru > 0 and eb_mask[ru - 1, cc]:
            ru -= 1
        while rd < H - 1 and eb_mask[rd + 1, cc]:
            rd += 1
        length = rd - ru + 1
        if length > best[0]:
            best = (length, cc, ru + r0, rd + r0)
    if best[0] == 0:
        return None
    # Refine stem column center at the midpoint row
    _, cc0, ru_abs, rd_abs = best
    rm = (ru_abs + rd_abs) // 2 - r0
    cl, cr = cc0, cc0
    while cl > 0 and eb_mask[rm, cl - 1]:
        cl -= 1
    while cr < eb_mask.shape[1] - 1 and eb_mask[rm, cr + 1]:
        cr += 1
    return 0.5 * (cl + cr) + 0.5, ru_abs, rd_abs


def read_fig6d(image: np.ndarray, variant: str, category: str) -> list[dict]:
    """Read one variant ('WT' or 'GABmut') and category ('shell', 'ring', 'disk') across Actin:VASP = 1:20, 1:10, 1:7.5."""
    info = pc.find_plot_box(image)
    y_axis = pc.fit_axis(info["y_ticks"], [100.0, 80.0, 60.0, 40.0, 20.0, 0.0])
    x_ticks = info["x_ticks"]
    if len(x_ticks) != 3:
        raise ValueError(f"Expected 3 x-ticks in Fig. 6D, got {len(x_ticks)}")
    c0, r0, c1, r1 = info["box"]
    px_per_pct = abs(y_axis.to_pixel(0.0) - y_axis.to_pixel(100.0)) / 100.0
    y_misfit = max(y_axis.residual_px, 0.5)
    err_pct = (2.0 + y_misfit) / px_per_pct
    eb_err_pct = (2.0 + y_misfit) / px_per_pct

    readings = []
    for ratio, xt in zip(RATIOS, x_ticks):
        if variant == "WT":
            a, b = int(round(xt - 182)), int(round(xt - 3))
        elif variant == "GABmut":
            a, b = int(round(xt + 3)), int(round(xt + 182))
        else:
            raise ValueError(f"Unknown variant {variant!r}")

        sub = image[r0:r1, a:b]
        is_red = (sub[:, :, 0] > 130) & (sub[:, :, 1] < 70) & (sub[:, :, 2] < 70)
        is_lg = (sub[:, :, 1] > 130) & (sub[:, :, 0] < 110) & (sub[:, :, 2] < 110)
        r_rows = np.flatnonzero(is_red.mean(axis=1) > 0.25) + r0
        lg_rows = np.flatnonzero(is_lg.mean(axis=1) > 0.25) + r0
        if not r_rows.size or not lg_rows.size:
            raise ValueError(f"Could not find bar segments for {variant} 1:{ratio:g}")
        row_shell = float(r_rows.min())
        row_ring = float(lg_rows.min())

        y_shell = min(100.0, max(0.0, y_axis.to_value(row_shell)))
        y_ring_cum = min(100.0, max(0.0, y_axis.to_value(row_ring)))
        y_ring = max(0.0, y_ring_cum - y_shell)
        y_disk = max(0.0, 100.0 - y_ring_cum)

        eb_mask = (sub.max(axis=2) < 55) & ((sub.max(axis=2) - sub.min(axis=2)) < 20)
        eb_s = _find_error_bar(eb_mask, r0, int(round(row_shell)), 72, 112)
        eb_r = _find_error_bar(eb_mask, r0, int(round(row_ring)), 40, 78)
        eb_d = _find_error_bar(eb_mask, r0, r0 + 8, 102, 148)

        if category == "shell":
            val = y_shell
            if eb_s is None:
                raise ValueError(f"Missing Shell error bar for {variant} 1:{ratio:g}")
            c_rel, ru, rd = eb_s
            yhi, ylo = y_axis.to_value(ru), y_axis.to_value(rd + 1)
            if ylo < 1.0:
                sd = yhi - y_shell
                clip_note = "lower end clipped at 0%; SD taken from upper half"
            elif yhi > 99.0:
                sd = y_shell - ylo
                clip_note = "upper end clipped at 100%; SD taken from lower half"
            else:
                sd = 0.5 * (yhi - ylo)
                clip_note = f"full symmetric SD bar [{ylo:.2f}, {yhi:.2f}]%"
            px_col, px_row = a + c_rel, row_shell
            note = f"Fig. 6D {variant} Shell at 1:{ratio:g}; SD={sd:.2f}% ({clip_note}; ry={y_axis.residual_px:.2f}px)"
        elif category == "ring":
            val = y_ring
            if eb_r is None:
                raise ValueError(f"Missing Ring error bar for {variant} 1:{ratio:g}")
            c_rel, ru, rd = eb_r
            yhi, ylo = y_axis.to_value(ru), y_axis.to_value(rd + 1)
            if yhi > 99.0:
                sd = y_ring_cum - ylo
                clip_note = "upper end clipped at 100%; SD taken from lower half"
            else:
                sd = 0.5 * (yhi - ylo)
                clip_note = f"full symmetric SD bar around cumulative top {y_ring_cum:.2f}%"
            px_col, px_row = a + c_rel, row_ring
            note = f"Fig. 6D {variant} Ring stacked height at 1:{ratio:g} (top={y_ring_cum:.2f}%); SD={sd:.2f}% ({clip_note}; ry={y_axis.residual_px:.2f}px)"
        elif category == "disk":
            val = y_disk
            if eb_d is None:
                raise ValueError(f"Missing Disc error bar for {variant} 1:{ratio:g}")
            c_rel, ru, rd = eb_d
            ylo = y_axis.to_value(rd + 1)
            sd = 100.0 - ylo
            px_col, px_row = a + c_rel, row_ring
            note = f"Fig. 6D {variant} Disc stacked height at 1:{ratio:g} (base={y_ring_cum:.2f}%); SD={sd:.2f}% (upper half clipped at 100%, lower end at {ylo:.2f}%; ry={y_axis.residual_px:.2f}px)"
        else:
            raise ValueError(f"Unknown category {category!r}")

        bar_lo = max(0.0, val - sd)
        bar_hi = min(100.0, val + sd)
        bar_lo_outer = max(0.0, bar_lo - eb_err_pct)
        bar_hi_outer = min(100.0, bar_hi + eb_err_pct)
        h_px = err_pct * px_per_pct
        readings.append({
            "x": float(ratio),
            "y": float(val),
            "lo": float(max(0.0, val - err_pct)),
            "hi": float(min(100.0, val + err_pct)),
            "bar_lo": float(bar_lo),
            "bar_hi": float(bar_hi),
            "bar_lo_outer": float(bar_lo_outer),
            "bar_hi_outer": float(bar_hi_outer),
            "px": [float(px_col), float(px_row)],
            "px_lo": float(px_row - h_px),
            "px_hi": float(px_row + h_px),
            "note": note,
        })
    return readings
