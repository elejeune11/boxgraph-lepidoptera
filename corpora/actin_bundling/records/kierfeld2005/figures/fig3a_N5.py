"""Calibrated reader for Fig. 3(a), N = 5 curve (open squares with dotted guide line)."""
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
    c_n5_drop = float(drops[2].mean())
    c_n3_drop = float(drops[3].mean())
    c_n2_drop = float(drops[4][0])

    out = []
    for x in at:
        col_f = mx * x + bx
        col = int(round(col_f))
        if col_f < c_n5_drop - 15:
            upper = _clusters_at(dark, col, r_top, r_zero - 20, hw=14, min_rows=10, min_mass=120)
            if not upper:
                out.append({"x": float(x), "y": None, "note": "no unbound plateau found"})
                continue
            rc, span, _, _ = upper[0]
            px_unc = max(10.0, 0.35 * span + 2.0 * y_res + 3.0)
        else:
            hw_use = 8 if x <= 5.3 else 16
            min_m = 80 if x <= 5.3 else 200
            lower = _clusters_at(dark, col, r_zero + 25, r_bot, hw=hw_use, min_rows=14, min_mass=min_m)
            if col_f < c_n3_drop - 15:
                idx = 0
            elif col_f < c_n2_drop - 15:
                idx = 1
            else:
                idx = 2
            # At large |W| >= 7.5, arrow 'b' can sit above N=5/N=10/N=20 planar bundle;
            # N=5 open squares sit at the bottom edge of the merged middle cluster (above N=10/20 compact).
            if x >= 7.4:
                # Restrict above arrow 'c' (below y = -12.8) so lower[-1] is always N=20 compact
                lower_clean = _clusters_at(dark, col, r_zero + 25, int(round(my * (-12.9) + by)), hw=16, min_rows=16, min_mass=200)
                # Middle bundle sits strictly above N=20 compact (lower_clean[-1]) and below N=3 (lower_clean[1])
                mid_candidates = [c for c in lower_clean if c[0] < lower_clean[-1][0] - 120 and c[0] > lower_clean[1][0] + 150]
                if not mid_candidates:
                    out.append({"x": float(x), "y": None, "note": "no middle bundle found"})
                    continue
                # Bottom-most of the middle bundle is the N=5 open-square edge
                _, span, g0, g1 = mid_candidates[-1]
                rc = g1 - 18.0 if span > 55 else 0.5 * (g0 + g1)
                px_unc = max(16.0, 2.0 * y_res + 14.0)
            else:
                if len(lower) <= idx:
                    out.append({"x": float(x), "y": None, "note": f"fewer than {idx+1} bound curves"})
                    continue
                rc, span, g0, g1 = lower[idx]
                if x >= 6.3 and span > 65:
                    # Merged briefly with planar branch around 6.5: take upper square center
                    rc = g0 + 22.0
                    px_unc = max(18.0, 2.0 * y_res + 16.0)
                else:
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
