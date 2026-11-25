"""Shared calibrated reader for Fig. 5—figure supplement 1 III(e) (PDF p. 42):
athermal relaxation curves S(t) vs t [s] for S_0 = 0.1 (blue), 0.2 (green),
0.3 (red), 0.4 (magenta), and their black dashed Eq. (12) exponential fits.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def _calibrate(image: np.ndarray):
    frame = C.find_box(image, dark=160, max_sat=80, min_frac=0.45)
    xt_c = C.detect_x_ticks(image, frame, side="out", spine="bottom", dark=160, max_sat=80, min_len=20, max_len=100)
    yt_c = C.detect_y_ticks(image, frame, side="out", spine="left", dark=160, max_sat=80, min_len=20, max_len=100)
    xt = C.select_regular_ticks(xt_c, expected_n=6, tol_px=15.0)
    yt = C.select_regular_ticks(yt_c, expected_n=5, tol_px=15.0)
    xa = C.calibrate_axis(xt, [0.0, 0.1, 0.2, 0.3, 0.4, 0.5], max_residual_px=4.0)
    ya = C.calibrate_axis(yt, [0.4, 0.3, 0.2, 0.1, 0.0], max_residual_px=4.0)
    return frame, xa, ya


def _color_mask(image: np.ndarray, s0: float) -> np.ndarray:
    R = image[:, :, 0].astype(int)
    G = image[:, :, 1].astype(int)
    B = image[:, :, 2].astype(int)
    if math.isclose(s0, 0.4):
        return (R > 140) & (B > 140) & (G < 70)
    if math.isclose(s0, 0.3):
        return (R > 170) & (G < 70) & (B < 70)
    if math.isclose(s0, 0.2):
        return (G > 75) & (R < 55) & (B < 55)
    if math.isclose(s0, 0.1):
        return (B > 170) & (R < 70) & (G < 70)
    raise ValueError(f"Unsupported S_0={s0}")


def _trace_sim_centerline(cmask: np.ndarray, box: tuple[float, float, float, float]) -> dict[int, float]:
    """Column-wise centerline of the colored simulation curve (used both for sim readings
    and as a clean narrow vertical guide band to isolate the co-located black dashed fit)."""
    c0, r0, c1, r1 = [int(v) for v in box]
    rows_by_col: dict[int, float] = {}
    for c in range(c0 + 8, c1 - 4):
        rr = np.where(cmask[r0 + 4 : r1 - 4, c])[0]
        if len(rr) >= 3:
            # Exclude very wide artifacts, take main run
            grps = C.groups(rr, gap=2)
            a, b = max(grps, key=lambda g: g[1] - g[0])
            rows_by_col[c] = (r0 + 4) + (a + b) / 2.0 + 0.5
    return rows_by_col


def read_sim(image: np.ndarray, at: list[float], s0: float) -> list[dict]:
    """Read solid colored simulation curve for target nematic order `s0`."""
    frame, xa, ya = _calibrate(image)
    cmask = _color_mask(image, s0)
    c0, r0, c1, r1 = [int(v) for v in frame["box"]]
    mask = np.zeros_like(cmask)
    mask[r0 + 4 : r1 - 4, c0 + 4 : c1 - 4] = cmask[r0 + 4 : r1 - 4, c0 + 4 : c1 - 4]

    # Where a black dash sits directly over the colored curve, interpolate smoothly across nearby
    # colored columns (window +-140 px, deg=2) so every requested x is read cleanly even under a black dash.
    return C.read_dashed_at(
        mask, xa, ya, at, frame["box"],
        window_px=150, sigma_px=45.0, deg=2, native_px=6.0,
        min_run_h=1.5, max_run_h=60.0,
        note_prefix=f"fig5s1e sim S0={s0:g}: ",
    )


def read_fit(image: np.ndarray, at: list[float], s0: float) -> list[dict]:
    """Read black dashed Eq. (12) exponential fit curve for target nematic order `s0`."""
    frame, xa, ya = _calibrate(image)
    cmask = _color_mask(image, s0)
    guide = _trace_sim_centerline(cmask, frame["box"])
    g_cols = np.array(sorted(guide.keys()), dtype=float)
    g_rows = np.array([guide[int(c)] for c in g_cols], dtype=float)
    if len(g_cols) < 20:
        raise RuntimeError(f"Fig 5-S1 III(e) S0={s0}: could not trace color guide")

    # Black dashed ink mask inside plot box
    R = image[:, :, 0].astype(int)
    G = image[:, :, 1].astype(int)
    B = image[:, :, 2].astype(int)
    blk = (R < 55) & (G < 55) & (B < 55)
    c0, r0, c1, r1 = [int(v) for v in frame["box"]]
    mask = np.zeros_like(blk)
    mask[r0 + 4 : r1 - 4, c0 + 4 : c1 - 4] = blk[r0 + 4 : r1 - 4, c0 + 4 : c1 - 4]

    # Restrict vertical search to +-38 px around the colored curve centerline
    # (curves are separated by ~480 px at plateau and >100 px for t >= 0.05 s; +-38 px cleanly excludes
    # both neighboring curves and the circled "(a)" annotation just below S_0 = 0.1 near t = 0.1 s).
    def band_fn(xv: float, col_px: float) -> tuple[float, float]:
        r_ref = float(np.interp(col_px, g_cols, g_rows))
        return (r_ref - 38.0, r_ref + 38.0)

    return C.read_dashed_at(
        mask, xa, ya, at, frame["box"],
        window_px=190, sigma_px=65.0, deg=2, native_px=6.0,
        band_fn=band_fn, min_run_h=2.0, max_run_h=50.0,
        note_prefix=f"fig5s1e fit S0={s0:g}: ",
    )
