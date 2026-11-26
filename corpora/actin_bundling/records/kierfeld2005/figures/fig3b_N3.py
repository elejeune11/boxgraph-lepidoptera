"""Calibrated reader for Fig. 3(b), N = 3 curve (pluses with dotted power-law guide line)."""
import numpy as np


def _groups(idx, gap=5):
    if len(idx) == 0:
        return []
    splits = np.where(np.diff(idx) > gap)[0] + 1
    return [g for g in np.split(idx, splits)]


def _calibrate(img):
    gray = img.astype(float).mean(axis=2)
    h, w = gray.shape
    dark = gray < 128
    vg = _groups(np.where(dark.sum(axis=0) > 0.5 * h)[0], 10)
    hg = _groups(np.where(dark.sum(axis=1) > 0.5 * w)[0], 10)
    x0, x1 = float(vg[0].mean()), float(vg[-1].mean())
    y0, y1 = float(hg[0].mean()), float(hg[-1].mean())

    top_t = _groups(np.where(dark[int(y0) + 12:int(y0) + 32, int(x0) + 50:int(x1) - 50].mean(axis=0) > 0.6)[0] + int(x0) + 50)
    x_px = np.array([x0, float(top_t[0].mean()), x1])
    x_log = np.array([-2.0, -1.0, 0.0])
    mx, bx = np.polyfit(x_log, x_px, 1)
    x_res = float(np.max(np.abs(x_px - (mx * x_log + bx))))

    rgt_t = _groups(np.where(dark[int(y0) + 50:int(y1) - 50, int(x1) - 32:int(x1) - 12].mean(axis=1) > 0.6)[0] + int(y0) + 50)
    y_px = np.array([y0, float(rgt_t[0].mean()), y1])
    y_log = np.array([0.0, -1.0, -2.0])
    my, by = np.polyfit(y_log, y_px, 1)
    y_res = float(np.max(np.abs(y_px - (my * y_log + by))))
    return (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res), dark


def _clusters_at(dark, col, r_lo, r_hi, hw=22, min_rows=8, min_mass=50):
    band = dark[r_lo:r_hi, col - hw:col + hw + 1]
    rr = np.where(band.sum(axis=1) >= 2)[0] + r_lo
    out = []
    for g in _groups(rr, 10):
        mass = int(band[g - r_lo, :].sum())
        if len(g) >= min_rows and mass >= min_mass:
            w = band[g - r_lo, :].sum(axis=1).astype(float)
            rc = float(np.sum(g * w) / np.sum(w))
            out.append((rc, float(g[-1] - g[0]), float(g[0]), float(g[-1]), mass))
    return out


def read(image, at):
    (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res), dark = _calibrate(image)
    r_bot = int(round(y1)) - 25

    out = []
    for x in at:
        col_f = mx * np.log10(x) + bx
        col = int(round(col_f))
        if x < 0.13:
            # Below N=5 (and below N=2 before their crossing near x = 0.14): N=3 is strictly the bottom-most marker
            # Exclude "b)" annotation near x in [0.015, 0.025], y < 0.04 (we only read x >= 0.037)
            r_lo = int(round(my * np.log10(0.098) + by))
            cl = _clusters_at(dark, col, r_lo, r_bot, hw=28, min_rows=8, min_mass=50)
            if not cl:
                out.append({"x": float(x), "y": None, "note": "no N=3 marker found"})
                continue
            rc, span, _, _, _ = cl[-1]
        elif x < 0.30:
            # After crossing N=2 around x = 0.14, N=3 plus markers lie strictly ABOVE N=2 diamonds
            # and slightly below the dotted guide line / text label "N=3"
            y_top = 0.028 if x < 0.22 else (0.021 if x < 0.26 else 0.0175)
            r_lo = int(round(my * np.log10(y_top) + by))
            cl = _clusters_at(dark, col, r_lo, r_bot, hw=22, min_rows=14, min_mass=110)
            if not cl:
                out.append({"x": float(x), "y": None, "note": "no N=3 marker found"})
                continue
            if len(cl) >= 2:
                rc, span, _, _, _ = cl[0]
            else:
                rc0, span, g0, g1, _ = cl[0]
                rc = g0 + 18.0 if span > 65 else rc0
        else:
            # For x >= 0.30, N=2 has ended; N=3 is strictly the bottom-most marker below the "N=3" text label
            r_lo = int(round(my * np.log10(0.018) + by))
            cl = _clusters_at(dark, col, r_lo, r_bot, hw=24, min_rows=8, min_mass=50)
            if not cl:
                out.append({"x": float(x), "y": None, "note": "no N=3 marker found"})
                continue
            rc, span, _, _, _ = cl[-1]
        px_unc = max(12.0, 0.35 * span + 2.0 * y_res + 4.0)
        y = 10.0 ** ((rc - by) / my)
        y_hi = 10.0 ** (((rc - px_unc) - by) / my)
        y_lo = 10.0 ** (((rc + px_unc) - by) / my)
        out.append({
            "x": float(x),
            "y": float(y),
            "lo": float(y_lo),
            "hi": float(y_hi),
            "px": [float(col_f), float(rc)],
            "px_lo": float(rc - px_unc),
            "px_hi": float(rc + px_unc),
            "note": f"x_res={x_res:.2f}px, y_res={y_res:.2f}px, span={span:.1f}px",
        })
    return out
