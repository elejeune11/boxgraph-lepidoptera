"""Calibrated reader for Fig. 2(b), blue dashed-triangle curve (Gamma_A = 13).
Uses 2D centroids of unoccluded blue dashes (plus bounding-box midpoint where a dash
straddles the black F=1 line) with symmetric local quadratic interpolation to avoid
triangle-centroid bias, beveled-dash sawtooth, and one-sided F=1 line occlusion."""
import numpy as np
from scipy.ndimage import label


def _groups(idx, gap=5):
    if len(idx) == 0:
        return []
    splits = np.where(np.diff(idx) > gap)[0] + 1
    return [g for g in np.split(idx, splits)]


def _calibrate(img):
    gray = img.astype(float).mean(axis=2)
    h, w = gray.shape
    vg = _groups(np.where((gray < 100).sum(axis=0) > 0.35 * h)[0], 10)
    hg = _groups(np.where((gray < 100).sum(axis=1) > 0.35 * w)[0], 10)
    x0, x1 = float(vg[0].mean()), float(vg[-1].mean())
    y0, y1 = float(hg[0].mean()), float(hg[-1].mean())

    top_t = _groups(np.where((gray[int(y0) + 6:int(y0) + 18, int(x0) + 25:int(x1) - 25] < 100).mean(axis=0) > 0.5)[0] + int(x0) + 25)
    x_px = np.array([x0, float(top_t[0].mean()), float(top_t[1].mean()), x1])
    x_val = np.array([0.0, 0.5, 1.5, 2.0])
    mx, bx = np.polyfit(x_val, x_px, 1)
    x_res = float(np.max(np.abs(x_px - (mx * x_val + bx))))

    lft_t = _groups(np.where((gray[int(y0) + 25:int(y1) - 25, int(x0) + 6:int(x0) + 18] < 100).mean(axis=1) > 0.5)[0] + int(y0) + 25)
    y_F1 = float(lft_t[3].mean())
    y_px = np.array([float(g.mean()) for g in lft_t] + [y1])
    y_val = np.array([6.0, 4.0, 2.0, 1.0, 0.0])
    my, by = np.polyfit(y_val, y_px, 1)
    y_res = float(np.max(np.abs(y_px - (my * y_val + by))))
    return (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res), y_F1


def read(image, at):
    (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res), y_F1 = _calibrate(image)
    blue = (image[:, :, 2] > 120) & (image[:, :, 0] < 80) & (image[:, :, 1] < 100)
    mask = np.zeros_like(blue)
    mask[int(y0) + 4:int(y1) - 4, int(x0) + 4:int(x1) - 4] = True
    mask[:int(y0 + 0.36 * (y1 - y0)), int(x0 + 0.55 * (x1 - x0)):] = False
    blue = blue & mask

    # Detect the F=0 cusp column where the curve touches the bottom border
    cusp_cols = np.nonzero(blue[int(y1) - 20:int(y1) - 4, int(x0) + 100:int(x0 + 0.45 * (x1 - x0))])[1]
    c_cusp = float(cusp_cols.mean()) + int(x0) + 100 if len(cusp_cols) > 0 else mx * 0.37 + bx

    # Extract 2D centroids of complete unoccluded blue dashes (area 380..550 px;
    # triangles have area ~850-950 px and clipped/occluded fragments have area < 350 px)
    lbl, n_comp = label(blue)
    dash_c, dash_r = [], []
    for i in range(1, n_comp + 1):
        rr, cc = np.where(lbl == i)
        if 380 <= len(rr) <= 550:
            dash_c.append(float(cc.mean()))
            dash_r.append(float(rr.mean()))

    # Reconstruct dashes cleanly straddling both sides of the horizontal black F=1 line
    str_c, str_r = [], []
    for c in range(int(c_cusp) + 20, int(x1) - 8):
        k = (c - bx) / mx
        if abs(k * 10.0 - round(k * 10.0)) < 0.22:
            continue
        rs = np.nonzero(blue[:, c])[0]
        if len(rs) >= 3 and rs.min() < y_F1 - 4 and rs.max() > y_F1 + 4 and 13 <= (rs.max() - rs.min() + 1) <= 18:
            str_c.append(float(c))
            str_r.append(0.5 * float(rs.min() + rs.max()))
    if str_c:
        str_arr = np.array(str_c)
        for g in _groups(np.array(str_c, dtype=int), gap=5):
            sel = np.isin(str_arr.astype(int), g)
            dash_c.append(float(str_arr[sel].mean()))
            dash_r.append(float(np.array(str_r)[sel].mean()))

    dash_c = np.array(dash_c)
    dash_r = np.array(dash_r)
    order = np.argsort(dash_c)
    dash_c, dash_r = dash_c[order], dash_r[order]

    out = []
    for x in at:
        col = mx * x + bx
        branch = (dash_c < c_cusp) if col < c_cusp else (dash_c > c_cusp)
        bc, br = dash_c[branch], dash_r[branch]
        # Pick 4 nearest dashes on this branch, ensuring col is bracketed when possible
        pick = np.argsort(np.abs(bc - col))[:min(4, len(bc))]
        if len(pick) == 4:
            if np.all(bc[pick] < col) and np.any(bc >= col):
                pick[-1] = int(np.where(bc >= col)[0][0])
            elif np.all(bc[pick] > col) and np.any(bc <= col):
                pick[-1] = int(np.where(bc <= col)[0][-1])
        if len(pick) < 2:
            out.append({"x": float(x), "y": None, "note": "insufficient clean blue dashes"})
            continue
        lc, lr = bc[pick], br[pick]
        deg = 2 if len(lc) >= 3 else 1
        poly = np.polyfit(lc - col, lr, deg)
        row = float(poly[-1])
        fit_res = float(np.sqrt(np.mean((np.polyval(poly, lc - col) - lr) ** 2)))
        y = (row - by) / my
        px_unc = max(5.0, 2.0 * y_res + 2.0 * fit_res + 3.0)
        dy = abs(px_unc / my)
        out.append({
            "x": float(x),
            "y": float(y),
            "lo": float(y - dy),
            "hi": float(y + dy),
            "px": [float(col), float(row)],
            "px_lo": float(row - px_unc),
            "px_hi": float(row + px_unc),
            "note": f"x_res={x_res:.2f}px, y_res={y_res:.2f}px, fit_res={fit_res:.2f}px",
        })
    return out
