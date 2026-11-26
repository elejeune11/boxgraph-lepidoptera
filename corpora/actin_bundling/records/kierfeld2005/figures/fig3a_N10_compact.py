"""Calibrated reader for Fig. 3(a), N = 10 compact lower branch (asterisks with dotted guide line)."""
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


def _clusters_at(dark, col, r_lo, r_hi, hw=8, min_rows=10, min_mass=60):
    band = dark[r_lo:r_hi, col - hw:col + hw + 1]
    rr = np.where(band.sum(axis=1) >= 2)[0] + r_lo
    out = []
    for g in _groups(rr, 8):
        mass = int(band[g - r_lo, :].sum())
        if len(g) >= min_rows and mass >= min_mass:
            w = band[g - r_lo, :].sum(axis=1).astype(float)
            rc = float(np.sum(g * w) / np.sum(w))
            out.append((rc, float(g[-1] - g[0]), float(g[0]), float(g[-1])))
    return out


def read(image, at):
    (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res), dark = _calibrate(image)
    r_zero = int(round(my * 0.0 + by))
    r_top = int(round(y0)) + 40
    r_bot = int(round(y1)) - 40

    r_d1 = int(round(my * (-0.3) + by))
    r_d2 = int(round(my * (-1.5) + by))
    col_sum = dark[r_d1:r_d2, int(x0) + 50:int(x0) + 1200].sum(axis=0)
    drops = _groups(np.where(col_sum > 12)[0] + int(x0) + 50, 15)
    c_n10_drop = float(drops[1].mean())

    out = []
    for x in at:
        col_f = mx * x + bx
        col = int(round(col_f))
        if col_f < c_n10_drop - 10:
            upper = _clusters_at(dark, col, r_top, r_zero - 20, hw=14, min_rows=10, min_mass=120)
            if not upper:
                out.append({"x": float(x), "y": None, "note": "no unbound plateau found"})
                continue
            rc, span, _, _ = upper[0]
            px_unc = max(10.0, 0.35 * span + 2.0 * y_res + 3.0)
        else:
            r_lim = int(round(my * (-9.0) + by)) if x < 5.1 else r_bot
            lower = _clusters_at(dark, col, r_zero + 25, r_lim, hw=8, min_rows=10, min_mass=60)
            if not lower:
                out.append({"x": float(x), "y": None, "note": "no N=10 compact marker found"})
                continue
            # Just after dropping (x in [4.72, 4.95]), N=10 compact is the topmost cluster below y=0
            # (above N=20 planar and N=20 compact). For x >= 5.1, it is the 2nd cluster from bottom
            # (immediately above N=20 compact, or upper half of their merged cluster near x=6.1).
            if x < 5.0:
                rc, span, _, _ = lower[0]
                px_unc = max(10.0, 0.35 * span + 2.0 * y_res + 3.0)
            else:
                last_rc, last_span, last_g0, last_g1 = lower[-1]
                if last_span > 65:
                    # Merged with N=20 compact diamond immediately below: asterisk sits in upper half
                    rc = last_g0 + 18.0
                    span = last_span
                    px_unc = max(14.0, 2.0 * y_res + 12.0)
                elif len(lower) >= 2:
                    rc, span, g0, g1 = lower[-2]
                    if span > 65:
                        # At x ~ 5.1, N=10 compact asterisk touches N=20 planar solid line directly above:
                        # take the lower asterisk center.
                        rc = g1 - 18.0
                        px_unc = max(14.0, 2.0 * y_res + 12.0)
                    else:
                        px_unc = max(10.0, 0.35 * span + 2.0 * y_res + 3.0)
                else:
                    out.append({"x": float(x), "y": None, "note": "N=10 compact not separated"})
                    continue
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
