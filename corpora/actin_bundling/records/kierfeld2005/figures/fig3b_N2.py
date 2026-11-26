"""Calibrated reader for Fig. 3(b), N = 2 curve (open diamonds with steepest dotted guide line)."""
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
    r_bot = int(round(y1)) - 20

    out = []
    for x in at:
        col_f = mx * np.log10(x) + bx
        col = int(round(col_f))
        if x < 0.025:
            # At x = 0.017, N=2 diamonds sit at the very top (~0.43), well above N=20 compact (~0.24)
            r_lo = int(round(y0)) + 40
            r_hi = int(round(my * np.log10(0.30) + by))
            cl = _clusters_at(dark, col, r_lo, r_hi, hw=28, min_rows=10, min_mass=70)
            if not cl:
                out.append({"x": float(x), "y": None, "note": "no top N=2 diamond found"})
                continue
            rc, span, _, _, _ = max(cl, key=lambda t: t[4])
        elif x < 0.06:
            # At x = 0.05, the N=2 diamond sits cleanly at ~0.11 between N=10 compact (~0.14) and N=5 (~0.095)
            r_lo = int(round(my * np.log10(0.125) + by))
            r_hi = int(round(my * np.log10(0.100) + by))
            cl = _clusters_at(dark, col, r_lo, r_hi, hw=24, min_rows=10, min_mass=70)
            if not cl:
                out.append({"x": float(x), "y": None, "note": "no N=2 diamond found near x=0.05"})
                continue
            rc, span, _, _, _ = max(cl, key=lambda t: t[4])
        elif x < 0.14:
            # On x in [0.08, 0.12], N=2 sits strictly below N=5 (~0.067..0.077) and above N=3 (~0.037..0.048)
            y_top = 0.068 if x < 0.095 else 0.056
            r_lo = int(round(my * np.log10(y_top) + by))
            r_hi = int(round(my * np.log10(0.032) + by))
            cl = _clusters_at(dark, col, r_lo, r_hi, hw=24, min_rows=10, min_mass=80)
            if not cl:
                out.append({"x": float(x), "y": None, "note": "no N=2 diamond found"})
                continue
            rc, span, _, _, _ = cl[0]
        else:
            # For x >= 0.18, N=2 has crossed below N=3 and is strictly the bottom-most marker curve
            # (above the "N=2" text label which sits around x in [0.18, 0.25], y < 0.0145)
            r_lo = int(round(my * np.log10(0.030) + by))
            y_floor = 0.0155 if x < 0.24 else 0.0105
            r_hi = int(round(my * np.log10(y_floor) + by))
            cl = _clusters_at(dark, col, r_lo, r_hi, hw=20, min_rows=10, min_mass=80)
            if not cl:
                out.append({"x": float(x), "y": None, "note": "no bottom N=2 diamond found"})
                continue
            rc0, span, g0, g1, _ = cl[-1]
            # If N=2 diamond touches N=3 plus right above, take lower diamond center
            rc = g1 - 18.0 if span > 65 else rc0
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
