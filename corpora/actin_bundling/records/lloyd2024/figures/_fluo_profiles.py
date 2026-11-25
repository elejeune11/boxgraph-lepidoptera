"""Shared calibrated reader for black-background TauSTED actin (green) & chitin (magenta) intensity profiles (Fig. 7E-F, SI Fig. 9F)."""

from __future__ import annotations

import numpy as np
from boxgraph import digitize as D
from _common import find_runs, run_centers, y_interval


def _trace_channel(image: np.ndarray, xa: D.Axis, ya: D.Axis,
                   r_top: int, r_bot: int, g_base: float,
                   channel: str, at: list[float], label: str) -> list[dict]:
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    H, W = R.shape

    green = (G > 140.0) & (G - R > 70.0) & (G - B > 70.0)
    magenta = (R > 130.0) & (B > 130.0) & (R - G > 60.0) & (B - G > 60.0)
    # Also detect white/bright colocalization pixels inside plot area if any (exclude text "(E)/(F)" at top-right)
    mask = green if channel == "actin" else magenta
    mask[:r_top, :] = False
    mask[r_bot:, :] = False
    green[:r_top, :] = False
    green[r_bot:, :] = False

    def col_row(m: np.ndarray, col: int) -> tuple[float | None, float]:
        if col < 0 or col >= W:
            return None, 0.0
        rs = np.where(m[:, col])[0]
        if len(rs) == 0:
            return None, 0.0
        # If there are multiple vertical segments due to a sharp spike, take the strongest/median run
        return float(np.median(rs) + 0.5), float(0.5 * (rs.max() - rs.min() + 1))

    out = []
    for xv in at:
        cf = float(xa.to_pixel(float(xv)))
        c = int(round(cf))
        rc, span = col_row(mask, c)
        note_extra = ""
        if rc is None:
            # Check small local neighborhood (+-3 px) first
            for dc in (1, -1, 2, -2, 3, -3):
                rc, span = col_row(mask, c + dc)
                if rc is not None:
                    break
        if rc is None:
            # Green is plotted on top of magenta in these TauSTED line profiles.
            # Look left and right up to 24 px for a narrow crossing; otherwise check baseline occlusion.
            rl, cl = None, None
            for dc in range(1, 26):
                r_try, _ = col_row(mask, c - dc)
                if r_try is not None:
                    rl, cl = r_try, c - dc
                    break
            rr, cr = None, None
            for dc in range(1, 26):
                r_try, _ = col_row(mask, c + dc)
                if r_try is not None:
                    rr, cr = r_try, c + dc
                    break
            if rl is not None and rr is not None and (cr - cl) <= 32:
                w = (cf - cl) / float(cr - cl)
                rc = (1.0 - w) * rl + w * rr
                span = 6.0
                note_extra = " [crossing interpolated under green]"
            else:
                # Occluded on baseline y = 0 under green curve
                rg, _ = col_row(green, c)
                if rg is not None and abs(rg - g_base) < 25.0:
                    rc = g_base
                    span = 4.0
                    note_extra = " [baseline occluded under green]"
                elif rl is not None and abs(rl - g_base) < 20.0:
                    rc = g_base
                    span = 4.0
                    note_extra = " [baseline]"
                elif rr is not None and abs(rr - g_base) < 20.0:
                    rc = g_base
                    span = 4.0
                    note_extra = " [baseline]"
                else:
                    raise RuntimeError(f"{label} ({channel}): no signal at x={xv} (col={c})")

        # Keep reading uncertainty proportional to line stroke / steepness (capped sensibly)
        half_px = max(3.5, min(span, 20.0))
        y_val, lo, hi, px_lo, px_hi = y_interval(ya, rc, half_px)
        # Clip tiny negative numerical jitter at baseline to >= 0 for non-negative intensity
        if abs(rc - g_base) < 2.0 and y_val < 0.0:
            y_val = 0.0
            lo = max(0.0, lo)
        out.append({
            "x": float(xv),
            "y": float(y_val), "lo": float(lo), "hi": float(hi),
            "px": [cf, float(rc)], "px_lo": float(px_lo), "px_hi": float(px_hi),
            "note": f"{label} {channel} (rx={xa.residual_px:.2f} px, ry={ya.residual_px:.2f} px){note_extra}",
        })
    return out


def read_fig7ef(image: np.ndarray, panel: str, channel: str, at: list[float]) -> list[dict]:
    """Read Fig. 7E (x=0..6 um, y=0..60) or Fig. 7F (x=0..8 um, y=0..60) from white-axis black-background plot."""
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    H, W = R.shape
    white = (R > 180.0) & (G > 180.0) & (B > 180.0)

    c0_band = np.where(white[180:H - 320, 100:420].mean(axis=0) > 0.70)[0] + 100
    c0_l, c0_r = int(c0_band.min()), int(c0_band.max())
    r1_band = np.where(white[H - 360:H - 50, c0_r + 20:W - 80].mean(axis=1) > 0.70)[0] + (H - 360)
    r1_t, r1_b = int(r1_band.min()), int(r1_band.max())

    # Green baseline row defines exact y = 0 level
    green = (G > 140.0) & (G - R > 70.0) & (G - B > 70.0)
    gb_rows = np.where(green[100:r1_t - 4, c0_r + 10:W - 60])[0] + 100
    g_base = float(np.median(gb_rows[gb_rows > gb_rows.max() - 15]) + 0.5)

    # Y digits (60, 40, 20, 0) sit strictly in [c0_l - 70 .. c0_l - 8], cleanly right of vertical 'Intensity' text
    y_runs = find_runs(np.where(white[100:r1_t, max(0, c0_l - 70):c0_l - 8].any(axis=1))[0] + 100, gap=12)
    if len(y_runs) != 4:
        raise RuntimeError(f"Fig. {panel}: expected 4 y-digit labels, found {len(y_runs)}: {y_runs}")
    y_cents = [r[2] for r in y_runs]
    # Align glyph-centroid baseline offset (3-5 px) cleanly to green y=0 baseline while preserving exact linear slope
    shift = g_base - y_cents[-1]
    if abs(shift) > 12.0:
        raise RuntimeError(f"Fig. {panel}: unexpected baseline shift {shift:.1f} px (y0={y_cents[-1]}, g_base={g_base})")
    yt = [yc + shift for yc in y_cents]

    # X digits sit strictly in [r1_b + 10 .. r1_b + 65], cleanly above 'Position (um)' label
    x_runs = find_runs(np.where(white[r1_b + 10:min(H, r1_b + 65), c0_r:W - 70].any(axis=0))[0] + c0_r, gap=12)
    x_lbls = [0.0, 2.0, 4.0, 6.0] if panel == "7E" else [0.0, 2.0, 4.0, 6.0, 8.0]
    if len(x_runs) != len(x_lbls):
        raise RuntimeError(f"Fig. {panel}: expected {len(x_lbls)} x-digit labels, found {len(x_runs)}: {x_runs}")
    xt = [r[2] for r in x_runs]

    xa = D.Axis("linear", list(zip(xt, x_lbls)))
    ya = D.Axis("linear", list(zip(yt, [60.0, 40.0, 20.0, 0.0])))
    return _trace_channel(image, xa, ya, r_top=40, r_bot=r1_t - 4, g_base=g_base,
                          channel=channel, at=at, label=f"Fig. {panel}")


def read_figs9f(image: np.ndarray, channel: str, at: list[float]) -> list[dict]:
    """Read SI Fig. 9F (x=0..10 um ticks, y=0..4000 ticks) from white-tick black-background plot."""
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    H, W = R.shape
    white = (R > 180.0) & (G > 180.0) & (B > 180.0)

    c0_band = np.where(white[120:H - 420, 350:680].mean(axis=0) > 0.70)[0] + 350
    c0_l, c0_r = int(c0_band.min()), int(c0_band.max())
    r1_band = np.where(white[H - 450:H - 250, c0_r + 20:W - 100].mean(axis=1) > 0.70)[0] + (H - 450)
    r1_t, r1_b = int(r1_band.min()), int(r1_band.max())

    yt = run_centers(np.where(white[80:r1_t, c0_l - 4])[0] + 80, gap=4)
    xt = run_centers(np.where(white[r1_b + 4, c0_r + 10:W - 200])[0] + c0_r + 10, gap=4)
    if len(yt) != 5 or len(xt) != 3:
        raise RuntimeError(f"SI Fig. 9F ticks mismatch: yt={yt}, xt={xt}")

    green = (G > 140.0) & (G - R > 70.0) & (G - B > 70.0)
    gb_rows = np.where(green[80:r1_t - 4, c0_r + 10:W - 80])[0] + 80
    g_base = float(np.median(gb_rows[gb_rows > gb_rows.max() - 15]) + 0.5)

    xa = D.Axis("linear", list(zip(xt, [0.0, 5.0, 10.0])))
    ya = D.Axis("linear", list(zip(yt, [4000.0, 3000.0, 2000.0, 1000.0, 0.0])))
    return _trace_channel(image, xa, ya, r_top=40, r_bot=r1_t - 4, g_base=g_base,
                          channel=channel, at=at, label="SI Fig. 9F")
