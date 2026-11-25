"""Calibrated reader for Fig. 3(a) Region A red '+' simulation markers phi_-(R) (p. 6, 1200 dpi)."""
from __future__ import annotations

import numpy as np
from scipy.ndimage import maximum_filter
from scipy.signal import fftconvolve


def _clusters(idxs: np.ndarray, gap: int = 3) -> list[tuple[float, int]]:
    if len(idxs) == 0:
        return []
    splits = np.where(np.diff(idxs) > gap)[0] + 1
    return [(float(g.mean()) + 0.5, int(len(g))) for g in np.split(idxs, splits)]


def _calibrate(img: np.ndarray):
    gray = img.astype(np.float32).mean(axis=2)
    dark = gray < 120.0
    r_lines = [c for c, _ in _clusters(np.where(dark.mean(axis=1) > 0.40)[0])]
    c_lines = [c for c, _ in _clusters(np.where(dark.mean(axis=0) > 0.40)[0])]
    r0, r1 = r_lines[0], r_lines[-1]
    c0, c1 = c_lines[0], c_lines[-1]
    ir0, ir1, ic0, ic1 = int(round(r0)), int(round(r1)), int(round(c0)), int(round(c1))

    bot_band = dark[ir1 - 18 : ir1 - 6, ic0 + 8 : ic1 - 8].mean(axis=0)
    xt_cands = [c + ic0 + 8 for c, w in _clusters(np.where(bot_band > 0.70)[0]) if w <= 12]
    x_val = np.array([0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0], dtype=float)
    x_px = np.array(
        [c0] + [min(xt_cands, key=lambda p: abs(p - (c0 + k * (c1 - c0) / 8.0))) for k in range(1, 8)] + [c1],
        dtype=float,
    )
    bx, ax = np.polyfit(x_px, x_val, 1)
    x_res = float(np.max(np.abs(x_px - (x_val - ax) / bx)))

    left_band = dark[ir0 + 8 : ir1 - 8, ic0 + 6 : ic0 + 18].mean(axis=1)
    yt_cands = [r + ir0 + 8 for r, w in _clusters(np.where(left_band > 0.70)[0]) if w <= 12]
    y_val = np.array([-0.7, -0.8, -0.9, -1.0], dtype=float)
    y_px = np.array(
        [r0] + [min(yt_cands, key=lambda p: abs(p - (r0 + k * (r1 - r0) / 3.0))) for k in range(1, 3)] + [r1],
        dtype=float,
    )
    by, ay = np.polyfit(y_px, y_val, 1)
    y_res = float(np.max(np.abs(y_px - (y_val - ay) / by)))
    return (ax, bx, x_res), (ay, by, y_res), (c0, r0, c1, r1)


def _plus_kernel(arm: int, thick: int) -> np.ndarray:
    s = 2 * arm + 1
    k = np.zeros((s, s), dtype=np.float32)
    c = arm
    k[:, c - thick // 2 : c + thick // 2 + 1] = 2.0
    k[c - thick // 2 : c + thick // 2 + 1, :] = 1.0
    for dr in (-1, 1):
        for dc in (-1, 1):
            r0, r1 = sorted([c + dr * (thick + 2), c + dr * (arm - 2)])
            c0, c1 = sorted([c + dc * (thick + 2), c + dc * (arm - 2)])
            k[r0:r1, c0:c1] = -0.3
    return k / k[k > 0].sum()


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    (ax, bx, x_res), (ay, by, y_res), (c0, r0, c1, r1) = _calibrate(image)
    R = image[:, :, 0].astype(np.int16)
    G = image[:, :, 1].astype(np.int16)
    B = image[:, :, 2].astype(np.int16)

    def x2c(v): return int(round((v - ax) / bx))
    def y2r(v): return int(round((v - ay) / by))

    mask = np.zeros(R.shape, dtype=bool)
    mask[int(r0) + 10 : int(r1) - 10, int(c0) + 10 : int(c1) - 10] = True
    mask[: y2r(-0.850), x2c(21.0) :] = False
    mask[y2r(-0.950) :, x2c(25.0) :] = False
    red = mask & (R - G > 70) & (R - B > 70)

    # 2x downsampled matched-filter detection for speed (~0.1s), followed by full-res refinement
    H, W = red.shape
    ds = red[: (H // 2) * 2, : (W // 2) * 2].reshape(H // 2, 2, W // 2, 2).max(axis=(1, 3)).astype(np.float32)
    kern_ds = _plus_kernel(arm=17, thick=4)
    score = fftconvolve(ds, kern_ds[::-1, ::-1], mode="same")
    mx = maximum_filter(score, size=21)
    prs, pcs = np.nonzero((score == mx) & (score >= 0.72))
    # Greedy NMS on score to eliminate flat-plateau ties
    cands = sorted([(float(score[r, c]), int(r) * 2 + 1, int(c) * 2 + 1) for r, c in zip(prs, pcs)], reverse=True)
    kept = []
    for s, r, c in cands:
        if all((r - kr) ** 2 + (c - kc) ** 2 > 36 ** 2 for _, kr, kc in kept):
            kept.append((s, r, c))

    kern_full = _plus_kernel(arm=34, thick=7)
    out = []
    for s_ds, r_est, c_est in sorted(kept, key=lambda t: t[2]):
        # Refine on full-resolution patch [r_est-46:r_est+47, c_est-46:c_est+47]
        pad = 46
        r_a, r_b = max(0, r_est - pad), min(H, r_est + pad + 1)
        c_a, c_b = max(0, c_est - pad), min(W, c_est + pad + 1)
        patch = red[r_a:r_b, c_a:c_b].astype(np.float32)
        sc = fftconvolve(patch, kern_full[::-1, ::-1], mode="same")
        # Restrict peak search to central +/- 10 px
        cr0, cc0 = r_est - r_a, c_est - c_a
        win = sc[max(0, cr0 - 10) : cr0 + 11, max(0, cc0 - 10) : cc0 + 11]
        dr_m, dc_m = np.unravel_index(int(np.argmax(win)), win.shape)
        rm = max(0, cr0 - 10) + dr_m
        cm = max(0, cc0 - 10) + dc_m
        sub = sc[max(0, rm - 3) : rm + 4, max(0, cm - 3) : cm + 4]
        w = np.maximum(0.0, sub - sub.min())
        yy, xx = np.mgrid[max(0, rm - 3) : rm + 4, max(0, cm - 3) : cm + 4]
        rc = r_a + float((yy * w).sum() / w.sum()) + 0.5
        cc = c_a + float((xx * w).sum() / w.sum()) + 0.5
        err_px = 4.0 + y_res
        x = float(ax + bx * cc)
        y = float(ay + by * rc)
        v1, v2 = float(ay + by * (rc - err_px)), float(ay + by * (rc + err_px))
        out.append({
            "x": x,
            "y": y,
            "lo": min(v1, v2),
            "hi": max(v1, v2),
            "px": [cc, rc],
            "px_lo": rc - err_px,
            "px_hi": rc + err_px,
            "note": f"red '+' marker (score={float(sc[rm, cm]):.2f}); cal residual ({x_res:.2f}, {y_res:.2f}) px",
        })
    return out
