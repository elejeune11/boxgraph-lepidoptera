"""Calibrated reader for Fig. 3(a), N = 3 curve (second bound branch from top, pluses)."""
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

    top_t = _groups(np.where(dark[int(y0) + 10:int(y0) + 30, int(x0) + 20:int(x1) - 20].mean(axis=0) > 0.6)[0] + int(x0) + 20)
    x_px = np.array([x0] + [float(g.mean()) for g in top_t] + [x1])
    x_val = np.linspace(4.0, 9.0, 11)
    mx, bx = np.polyfit(x_val, x_px, 1)
    x_res = float(np.max(np.abs(x_px - (mx * x_val + bx))))

    rgt_t = _groups(np.where(dark[int(y0) + 20:int(y1) - 20, int(x1) - 30:int(x1) - 10].mean(axis=1) > 0.6)[0] + int(y0) + 20)
    y_px = np.array([y0] + [float(g.mean()) for g in rgt_t] + [y1])
    y_val = np.array([2.0, 0.0, -2.0, -4.0, -6.0, -8.0, -10.0, -12.0, -14.0])
    my, by = np.polyfit(y_val, y_px, 1)
    y_res = float(np.max(np.abs(y_px - (my * y_val + by))))
    return (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res), dark


def _clusters_at(dark, col, r_lo, r_hi, hw=16, min_rows=14, min_mass=180):
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
    r_zero = int(round(my * 0.0 + by))
    r_top = int(round(y0)) + 40
    r_bot = int(round(y1)) - 40
    # Locate the N=3 vertical drop line directly from the image (4th vertical dotted line from left)
    r_d1 = int(round(my * (-0.3) + by))
    r_d2 = int(round(my * (-1.5) + by))
    col_sum = dark[r_d1:r_d2, int(x0) + 50:int(x0) + 1200].sum(axis=0)
    drops = _groups(np.where(col_sum > 12)[0] + int(x0) + 50, 15)
    # drops[0]=N20, drops[1]=N10, drops[2]=N5, drops[3]=N3, drops[4]=N2
    c_n3_drop = float(drops[3].mean())
    c_n2_drop = float(drops[4][0])

    out = []
    for x in at:
        col_f = mx * x + bx
        col = int(round(col_f))
        if col_f < c_n3_drop - 20:
            upper = _clusters_at(dark, col, r_top, r_zero - 20, hw=14, min_rows=10, min_mass=120)
            if not upper:
                out.append({"x": float(x), "y": None, "note": "no unbound plateau found"})
                continue
            rc, span, _ = upper[0]
        else:
            lower = _clusters_at(dark, col, r_zero + 25, r_bot, hw=18, min_rows=14, min_mass=150)
            # Before N=2 drops, N=3 is the topmost bound curve; after N=2 drops, N=3 is the 2nd bound curve!
            idx = 0 if col_f < c_n2_drop - 15 else 1
            if len(lower) <= idx:
                out.append({"x": float(x), "y": None, "note": f"fewer than {idx+1} bound curves"})
                continue
            rc, span, _ = lower[idx]
        px_unc = max(10.0, 0.35 * span + 2.0 * y_res + 3.0)
        y = (rc - by) / my
        dy = abs(px_unc / my)
        out.append({
            "x": float(x),
            "y": float(y),
            "lo": float(y - dy),
            "hi": float(y + dy),
            "px": [float(col_f), float(rc)],
            "px_lo": float(rc - px_unc),
            "px_hi": float(rc + px_unc),
            "note": f"x_res={x_res:.2f}px, y_res={y_res:.2f}px, span={span:.1f}px",
        })
    return out
