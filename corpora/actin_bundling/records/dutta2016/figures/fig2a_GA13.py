"""Calibrated reader for Fig. 2(a), blue dashed-triangle curve (Gamma_A = 13).
Samples on pure-line dashes (avoiding upward-triangle centroid bias) using local linear fit."""
import numpy as np


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
    x_px = np.array([x0] + [float(g.mean()) for g in top_t])
    x_val = np.array([0.0, 0.5, 1.0, 1.5])
    mx, bx = np.polyfit(x_val, x_px, 1)
    x_res = float(np.max(np.abs(x_px - (mx * x_val + bx))))

    rgt_t = _groups(np.where((gray[int(y0) + 25:int(y1) - 25, int(x1) - 18:int(x1) - 6] < 100).mean(axis=1) > 0.5)[0] + int(y0) + 25)
    y_top_tick = float(rgt_t[0].mean())
    y_bot_tick = float(rgt_t[-1].mean())
    mid_strip = gray[int(y0) + 200:int(y1) - 200, int(x0 + 0.3 * (x1 - x0)):int(x0 + 0.9 * (x1 - x0))]
    y_zero = float((np.where((mid_strip < 100).mean(axis=1) > 0.8)[0] + int(y0) + 200).mean())
    y_px = np.array([y_top_tick, y_zero, y_bot_tick])
    y_val = np.array([4.0, 0.0, -4.0])
    my, by = np.polyfit(y_val, y_px, 1)
    y_res = float(np.max(np.abs(y_px - (my * y_val + by))))
    return (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res)


def read(image, at):
    (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res) = _calibrate(image)
    blue = (image[:, :, 2] > 120) & (image[:, :, 0] < 80) & (image[:, :, 1] < 100)
    mask = np.zeros_like(blue)
    mask[int(y0) + 4:int(y1) - 4, int(x0) + 4:int(x1) - 4] = True
    mask[:int(y0 + 0.36 * (y1 - y0)), int(x0 + 0.55 * (x1 - x0)):] = False
    blue = blue & mask

    # Collect pure-line blue columns (excluding triangle markers at 0.1 multiples)
    line_c, line_r = [], []
    for c in range(int(x0) + 8, int(x1) - 8):
        k = (c - bx) / mx
        if abs(k * 10.0 - round(k * 10.0)) < 0.22:
            continue
        rows = np.nonzero(blue[:, c])[0]
        if len(rows) > 0 and (rows.max() - rows.min()) < 22:
            line_c.append(float(c))
            line_r.append(float(rows.mean()))
    line_c = np.array(line_c)
    line_r = np.array(line_r)

    out = []
    for x in at:
        col = mx * x + bx
        sel = np.abs(line_c - col) <= 45.0
        if sel.sum() < 4:
            sel = np.abs(line_c - col) <= 90.0
        if sel.sum() < 2:
            out.append({"x": float(x), "y": None, "note": "insufficient blue dash pixels"})
            continue
        lc, lr = line_c[sel], line_r[sel]
        deg = 2 if (sel.sum() >= 12 and (lc.max() - lc.min()) > 30) else 1
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
