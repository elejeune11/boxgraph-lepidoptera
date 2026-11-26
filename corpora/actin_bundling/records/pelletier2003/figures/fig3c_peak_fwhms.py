"""Calibrated reader for Fig. 3(c) (p. 3): full widths at half maximum (FWHM_1, FWHM_2) of plotted Gaussians G_1 and G_2."""

from __future__ import annotations

import numpy as np
from scipy.optimize import curve_fit

from boxgraph.digitize import Axis


def _cluster(idxs: np.ndarray, gap: int = 4) -> list[float]:
    if len(idxs) == 0:
        return []
    groups = np.split(idxs, np.where(np.diff(idxs) > gap)[0] + 1)
    return [float(g.mean()) for g in groups if len(g) > 0]


def _calibrate_fig3c(gray: np.ndarray):
    dark = gray < 140
    r1 = _cluster(np.where(dark[100:1700, 200:950].sum(axis=1) > 500)[0])[-1] + 100.0
    c0 = _cluster(np.where(dark[200:1550, 100:900].sum(axis=0) > 800)[0])[0] + 100.0

    xt = [x + int(c0) - 15 for x in _cluster(np.where(gray[int(r1) + 8 : int(r1) + 18, int(c0) - 15 : 950].mean(axis=0) < 200)[0])]
    ax_x = Axis("linear", list(zip(xt[:4], [0.1, 0.2, 0.3, 0.4])))

    yt = [y + 200.0 for y in _cluster(np.where(gray[200 : int(r1) + 10, int(c0) - 16 : int(c0) - 6].mean(axis=1) < 200)[0])]
    ax_y = Axis("linear", list(zip(yt[:5], [800.0, 600.0, 400.0, 200.0, 0.0])))
    return ax_x, ax_y, (c0, r1)


def _trace_gaussians(gray: np.ndarray, ax_x: Axis, ax_y: Axis, r1: float):
    g1_q, g1_y = [], []
    for c in range(int(ax_x.to_pixel(0.181)), int(ax_x.to_pixel(0.203)), 4):
        rr = np.where(gray[540:760, c] < 160)[0] + 540
        if len(rr):
            seg = np.split(rr, np.where(np.diff(rr) > 5)[0] + 1)[-1]
            g1_q.append(ax_x.to_value(c))
            g1_y.append(ax_y.to_value(float(seg.mean())))

    for yv in np.arange(580, 460, -20):
        r = int(round(ax_y.to_pixel(yv)))
        cc = np.where(gray[r, 370:480] < 160)[0] + 370
        segs = np.split(cc, np.where(np.diff(cc) > 5)[0] + 1) if len(cc) else []
        thin = [s for s in segs if len(s) <= 16]
        lefts = [s for s in thin if 0.165 < ax_x.to_value(float(s.mean())) < 0.190]
        rights = [s for s in thin if 0.192 < ax_x.to_value(float(s.mean())) < 0.220]
        if lefts:
            s = max(lefts, key=lambda x: x.mean())
            g1_q.append(ax_x.to_value(float(s.mean())))
            g1_y.append(float(yv))
        if rights:
            s = min(rights, key=lambda x: x.mean())
            g1_q.append(ax_x.to_value(float(s.mean())))
            g1_y.append(float(yv))

    for yv in np.arange(440, 60, -20):
        if 240 <= yv <= 290:
            continue
        r = int(round(ax_y.to_pixel(yv)))
        cc = np.where(gray[r, 455:560] < 160)[0] + 455
        segs = np.split(cc, np.where(np.diff(cc) > 5)[0] + 1) if len(cc) else []
        thin = [s for s in segs if len(s) <= 14]
        if yv > 290 and thin:
            g1_q.append(ax_x.to_value(float(thin[0].mean())))
            g1_y.append(float(yv))
        elif yv < 240 and thin:
            g1_q.append(ax_x.to_value(float(thin[-1].mean())))
            g1_y.append(float(yv))

    g2_q, g2_y = [], []
    for yv in np.arange(60, 235, 15):
        r = int(round(ax_y.to_pixel(yv)))
        cc = np.where(gray[r, 340:500] < 160)[0] + 340
        segs = np.split(cc, np.where(np.diff(cc) > 5)[0] + 1) if len(cc) else []
        cands = [s for s in segs if len(s) <= 14 and 0.160 < ax_x.to_value(float(s.mean())) < 0.222]
        if cands:
            g2_q.append(ax_x.to_value(float(cands[0].mean())))
            g2_y.append(float(yv))

    for yv in np.arange(300, 415, 15):
        r = int(round(ax_y.to_pixel(yv)))
        cc = np.where(gray[r, 500:600] < 160)[0] + 500
        segs = np.split(cc, np.where(np.diff(cc) > 5)[0] + 1) if len(cc) else []
        cands = [s for s in segs if len(s) <= 16]
        if cands:
            g2_q.append(ax_x.to_value(float(cands[-1].mean())))
            g2_y.append(float(yv))

    for q in np.arange(0.285, 0.455, 0.01):
        c = int(round(ax_x.to_pixel(q)))
        rr = np.where(gray[560 : int(r1) - 8, c] < 150)[0] + 560
        if len(rr):
            seg = np.split(rr, np.where(np.diff(rr) > 5)[0] + 1)[0]
            g2_q.append(float(q))
            g2_y.append(ax_y.to_value(float(seg.mean())))

    def _one_g(q_arr, A, mu, sigma):
        return A * np.exp(-0.5 * ((q_arr - mu) / sigma) ** 2)

    p1, cov1 = curve_fit(_one_g, np.array(g1_q), np.array(g1_y), p0=[613.0, 0.19, 0.028])
    p2, cov2 = curve_fit(_one_g, np.array(g2_q), np.array(g2_y), p0=[445.0, 0.29, 0.066])
    return (p1, np.sqrt(np.diag(cov1))), (p2, np.sqrt(np.diag(cov2)))


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    gray = image.astype(float).mean(axis=2)
    ax_x, ax_y, (_, r1) = _calibrate_fig3c(gray)
    (p1, se1), (p2, se2) = _trace_gaussians(gray, ax_x, ax_y, r1)

    k_fwhm = 2.0 * np.sqrt(2.0 * np.log(2.0))
    px_quant = 2.0
    # Two-sided width combines calibration scale error + line thickness (~2px) + fit uncertainty
    cal_width_err = float(np.sqrt(2.0) * np.hypot(px_quant, ax_x.residual_px) * abs(ax_x.b))

    entries = [
        (p1[0], p1[1], abs(p1[2]) * k_fwhm, se1[2] * k_fwhm, 2.0 * abs(p1[2]), "G_1 FWHM"),
        (p2[0], p2[1], abs(p2[2]) * k_fwhm, se2[2] * k_fwhm, 2.0 * abs(p2[2]), "G_2 FWHM"),
    ]
    out = []
    for q_at in (at or [0.19, 0.29]):
        # Match to nearest Gaussian peak (0.19 -> G_1, 0.29 -> G_2)
        A, mu, fwhm, se_fwhm, w2s, desc = min(entries, key=lambda e: abs(e[1] - q_at))
        df = float(np.hypot(cal_width_err, 1.96 * se_fwhm))
        # Mark right half-maximum point on the plotted Gaussian curve for visual verification
        c_hm = float(ax_x.to_pixel(mu + 0.5 * fwhm))
        r_hm = float(ax_y.to_pixel(0.5 * A))
        out.append({
            "x": float(q_at),
            "y": float(fwhm),
            "lo": float(fwhm - df),
            "hi": float(fwhm + df),
            "px": [c_hm, r_hm],
            "note": f"{desc}={fwhm:.4f} (2*sigma={w2s:.4f}) nm^-1; x_res={ax_x.residual_px:.2f}px",
        })
    return out
