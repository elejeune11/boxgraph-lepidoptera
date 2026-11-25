"""Shared calibrated reader for Fig. 5B, 5D, and 5E (8-curve CytoSim plots across [Arp2/3]=0..150 nM)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C

ARP_CONCS = [0, 10, 25, 50, 75, 100, 125, 150]
ARP_RGB = [
    np.array([85.0, 149.0, 78.0], dtype=np.float32),    # 0 nM
    np.array([116.0, 163.0, 101.0], dtype=np.float32),  # 10 nM
    np.array([168.0, 205.0, 140.0], dtype=np.float32),  # 25 nM
    np.array([194.0, 159.0, 186.0], dtype=np.float32),  # 50 nM
    np.array([173.0, 122.0, 160.0], dtype=np.float32),  # 75 nM
    np.array([154.0, 84.0, 137.0], dtype=np.float32),   # 100 nM
    np.array([138.0, 39.0, 119.0], dtype=np.float32),   # 125 nM
    np.array([95.0, 34.0, 86.0], dtype=np.float32),     # 150 nM
]


def _color_mask_for_conc(image: np.ndarray, conc_nM: int, r0: int, r1: int, c0: int, c1: int, tol: float = 18.0) -> np.ndarray:
    """Fast nearest-palette mask inside [r0:r1, c0:c1] matching only solid-line pixels of `conc_nM`."""
    idx = ARP_CONCS.index(conc_nM)
    sub = image[r0:r1, c0:c1].astype(np.float32)
    R, G, B = sub[:, :, 0], sub[:, :, 1], sub[:, :, 2]
    # Compute squared distance to target color first, then check it is strictly closest among all 8 palette colors
    t = ARP_RGB[idx]
    d2_t = (R - t[0]) ** 2 + (G - t[1]) ** 2 + (B - t[2]) ** 2
    cand = d2_t < (tol ** 2)
    for j, col in enumerate(ARP_RGB):
        if j == idx:
            continue
        d2_j = (R - col[0]) ** 2 + (G - col[1]) ** 2 + (B - col[2]) ** 2
        cand &= (d2_t <= d2_j)
    full = np.zeros(image.shape[:2], dtype=bool)
    full[r0:r1, c0:c1] = cand
    return full


def read_fig5b(image: np.ndarray, at: list[float], conc_nM: int) -> list[dict]:
    dark, gray, hl, vl = C.find_l_spines(image, thresh=110.0, row_frac=0.35, col_frac=0.35)
    v0 = vl[0]
    y_end = np.where(dark[:, int(round(v0[0]))])[0].max()
    hb = min(hl, key=lambda h: abs(h[0] - y_end))

    yt_c = C.find_outward_ticks(dark, v0, (0, int(hb[2]) + 25), "y", side=-1, band=(5, 35), frac=0.45, max_width=24)
    yt = C.match_ticks_regular(yt_c, 7, tol=8.0)
    ya = C.verify_axis("linear", list(zip(yt, [6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0])), max_residual_px=2.5)

    xt_c = C.find_outward_ticks(dark, hb, (int(v0[1]) - 25, image.shape[1] - 5), "x", side=+1, band=(5, 35), frac=0.45, max_width=24)
    # 5 regular outward ticks correspond to 0.2, 0.4, 0.6, 0.8, 1.0 um (xt_c[0] is the left spine near 0.11)
    xt = C.match_ticks_regular(xt_c[1:] if len(xt_c) >= 6 else xt_c, 5, tol=8.0)
    xa = C.verify_axis("linear", list(zip(xt, [0.2, 0.4, 0.6, 0.8, 1.0])), max_residual_px=2.5)

    r_min, r_max = int(yt[0]) - 20, int(hb[2]) + 10
    c_min, c_max = int(v0[2]) + 2, int(xt[-1]) + 25
    mask = _color_mask_for_conc(image, conc_nM, r_min, r_max, c_min, c_max, tol=18.0)
    return [
        C.sample_solid_curve_y(mask, xa, ya, float(x), r_min, r_max, search_dx=5, note_prefix=f"Fig5B [Arp2/3]={conc_nM}nM, ")
        for x in at
    ]


def read_fig5d(image: np.ndarray, at: list[float], conc_nM: int) -> list[dict]:
    dark, gray, hl, vl = C.find_l_spines(image, thresh=110.0, row_frac=0.35, col_frac=0.35)
    v0 = vl[0]
    y_end = np.where(dark[:, int(round(v0[0]))])[0].max()
    hb = min(hl, key=lambda h: abs(h[0] - y_end))

    yt_c = C.find_outward_ticks(dark, v0, (0, int(hb[2]) + 25), "y", side=-1, band=(5, 35), frac=0.45, max_width=24)
    yt = C.match_ticks_regular(yt_c, 6, tol=8.0)
    ya = C.verify_axis("linear", list(zip(yt, [1000.0, 800.0, 600.0, 400.0, 200.0, 0.0])), max_residual_px=2.5)

    xt_c = C.find_outward_ticks(dark, hb, (int(v0[1]) - 25, image.shape[1] - 5), "x", side=+1, band=(5, 35), frac=0.45, max_width=24)
    xt = C.match_ticks_regular(xt_c, 4, tol=8.0)
    xa = C.verify_axis("linear", list(zip(xt, [0.0, 200.0, 400.0, 600.0])), max_residual_px=2.5)

    r_min, r_max = max(0, int(yt[0]) - 60), int(hb[1]) - 2
    c_min, c_max = int(v0[2]) + 2, int(xt[-1]) + 25
    # Tight tolerance (17) isolates the opaque solid mean curve and excludes pale transparent SD shading
    mask = _color_mask_for_conc(image, conc_nM, r_min, r_max, c_min, c_max, tol=17.0)
    return [
        C.sample_solid_curve_y(mask, xa, ya, float(x), r_min, r_max, search_dx=5, note_prefix=f"Fig5D [Arp2/3]={conc_nM}nM, ")
        for x in at
    ]


def read_fig5e(image: np.ndarray, at: list[float], conc_nM: int) -> list[dict]:
    dark, gray, hl, vl = C.find_l_spines(image, thresh=110.0, row_frac=0.35, col_frac=0.35)
    v0 = vl[0]
    y_end = np.where(dark[:, int(round(v0[0]))])[0].max()
    hb = min(hl, key=lambda h: abs(h[0] - y_end))

    yt_c = C.find_outward_ticks(dark, v0, (0, int(hb[2]) + 25), "y", side=-1, band=(5, 35), frac=0.45, max_width=24)
    yt = C.match_ticks_regular(yt_c, 5, tol=8.0)
    ya = C.verify_axis("linear", list(zip(yt, [1.0, 0.75, 0.50, 0.25, 0.0])), max_residual_px=2.5)

    xt_c = C.find_outward_ticks(dark, hb, (int(v0[1]) - 25, image.shape[1] - 5), "x", side=+1, band=(5, 35), frac=0.45, max_width=24)
    xt = C.match_ticks_regular(xt_c, 4, tol=8.0)
    xa = C.verify_axis("linear", list(zip(xt, [0.0, 2.0, 4.0, 6.0])), max_residual_px=2.5)

    r_min, r_max = max(0, int(yt[0]) - 25), int(hb[2]) + 10
    c_min, c_max = int(v0[2]) + 2, image.shape[1] - 5
    mask = _color_mask_for_conc(image, conc_nM, r_min, r_max, c_min, c_max, tol=18.0)
    return [
        C.sample_solid_curve_y(mask, xa, ya, float(x), r_min, r_max, search_dx=5, note_prefix=f"Fig5E [Arp2/3]={conc_nM}nM, ")
        for x in at
    ]
