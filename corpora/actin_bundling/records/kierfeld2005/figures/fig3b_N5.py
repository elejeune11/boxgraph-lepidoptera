"""Calibrated reader for Fig. 3(b), N = 5 curve (open squares with dotted guide line)."""
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


def _clusters_at(dark, col, r_lo, r_hi, hw=22, min_rows=10, min_mass=80):
    band = dark[r_lo:r_hi, col - hw:col + hw + 1]
    rr = np.where(band.sum(axis=1) >= 2)[0] + r_lo
    out = []
    for g in _groups(rr, 12):
        mass = int(band[g - r_lo, :].sum())
        if len(g) >= min_rows and mass >= min_mass:
            w = band[g - r_lo, :].sum(axis=1).astype(float)
            rc = float(np.sum(g * w) / np.sum(w))
            out.append((rc, float(g[-1] - g[0]), mass))
    return out


def read(image, at):
    (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res), dark = _calibrate(image)
    out = []
    for x in at:
        col_f = mx * np.log10(x) + bx
        col = int(round(col_f))
        # Open-square N=5 band lies cleanly between N=10 compact above and N=2/N=3 below
        if x < 0.05:
            y_max, y_min = 0.155, 0.098
        elif x < 0.13:
            y_max, y_min = 0.095, 0.062
        else:
            y_max, y_min = 0.070, 0.042
        r_lo = int(round(my * np.log10(y_max) + by))
        r_hi = int(round(my * np.log10(y_min) + by))
        cl = _clusters_at(dark, col, r_lo, r_hi, hw=22, min_rows=14, min_mass=180)
        if not cl:
            out.append({"x": float(x), "y": None, "note": "no N=5 square marker found"})
            continue
        # At x = 0.04, N=10 compact sits at ~0.15 above N=5 (~0.108), so take the bottom-most square cluster.
        # At x = 0.10, two stacked squares sit around 0.071..0.081; take lower square on the guide line.
        rc, span, _ = cl[-1] if (x < 0.05 or x == 0.10) else max(cl, key=lambda t: t[2])
        px_unc = max(10.0, 0.35 * span + 2.0 * y_res + 4.0)
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
