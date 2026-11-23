"""Calibrated reader for Fig. 2(a), black solid-square curve (Gamma_A = 31)."""
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
    return (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res), y_zero


def read(image, at):
    (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res), y_zero = _calibrate(image)
    # Neutral dark black pixels strictly below the y=0 line (lowest curve in Fig. 2(a))
    r_ch = image[:, :, 0].astype(float)
    g_ch = image[:, :, 1].astype(float)
    b_ch = image[:, :, 2].astype(float)
    blk = (r_ch < 60) & (g_ch < 60) & (b_ch < 60)

    out = []
    for x in at:
        col = int(round(mx * x + bx))
        strip = blk[int(y_zero) + 12:int(y1) - 35, max(int(x0) + 4, col - 3):min(int(x1) - 4, col + 4)]
        rows = np.where(strip.mean(axis=1) > 0.5)[0] + int(y_zero) + 12
        if len(rows) == 0:
            out.append({"x": float(x), "y": None, "note": "no black curve pixels below y=0"})
            continue
        # Lowest cluster is the Gamma_A = 31 curve (squares are vertically symmetric)
        grps = _groups(rows, gap=6)
        g = grps[-1]
        row = float(g.mean())
        y = (row - by) / my
        px_unc = max(4.5, 2.0 * y_res + 3.0)
        dy = abs(px_unc / my)
        out.append({
            "x": float(x),
            "y": float(y),
            "lo": float(y - dy),
            "hi": float(y + dy),
            "px": [float(col), float(row)],
            "px_lo": float(row - px_unc),
            "px_hi": float(row + px_unc),
            "note": f"x_res={x_res:.2f}px, y_res={y_res:.2f}px",
        })
    return out
