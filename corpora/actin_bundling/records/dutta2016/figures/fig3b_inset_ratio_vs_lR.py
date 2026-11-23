"""Calibrated reader for Fig. 3(b) inset, solid red curve (Bundle size / lambda_R vs. lambda_R)."""
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

    bot_t = _groups(np.where((gray[int(y1) - 18:int(y1) - 6, int(x0) + 25:int(x1) - 25] < 100).mean(axis=0) > 0.5)[0] + int(x0) + 25)
    if len(bot_t) != 3:
        raise RuntimeError(f"Expected 3 x-ticks (0.5,0.7,0.9) in Fig. 3(b) inset, found {len(bot_t)}")
    x_px = np.array([float(g.mean()) for g in bot_t])
    x_val = np.array([0.5, 0.7, 0.9])
    mx, bx = np.polyfit(x_val, x_px, 1)
    x_res = float(np.max(np.abs(x_px - (mx * x_val + bx))))

    lft_t = _groups(np.where((gray[int(y0) + 25:int(y1) - 25, int(x0) + 6:int(x0) + 18] < 100).mean(axis=1) > 0.5)[0] + int(y0) + 25)
    if len(lft_t) != 3:
        raise RuntimeError(f"Expected 3 y-ticks (30,20,10) in Fig. 3(b) inset, found {len(lft_t)}")
    y_px = np.array([float(g.mean()) for g in lft_t])
    y_val = np.array([30.0, 20.0, 10.0])
    my, by = np.polyfit(y_val, y_px, 1)
    y_res = float(np.max(np.abs(y_px - (my * y_val + by))))
    return (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res)


def read(image, at):
    (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res) = _calibrate(image)
    red = (image[:, :, 0] > 180) & (image[:, :, 1] < 80) & (image[:, :, 2] < 80)
    mask = np.zeros_like(red)
    mask[int(y0) + 4:int(y1) - 4, int(x0) + 4:int(x1) - 4] = True
    red = red & mask

    out = []
    for x in at:
        col = int(round(mx * x + bx))
        rows = np.nonzero(red[:, col])[0]
        if len(rows) == 0:
            win = red[:, max(int(x0) + 4, col - 3):min(int(x1) - 4, col + 4)]
            rows, _ = np.nonzero(win)
        if len(rows) == 0:
            out.append({"x": float(x), "y": None, "note": "no red curve pixels found"})
            continue
        row = float(rows.mean())
        half_thick = 0.5 * float(rows.max() - rows.min())
        slope_px = 0.0
        r_prev = np.nonzero(red[:, max(int(x0) + 4, col - 5)])[0]
        r_next = np.nonzero(red[:, min(int(x1) - 4, col + 5)])[0]
        if len(r_prev) > 0 and len(r_next) > 0:
            slope_px = abs(float(r_next.mean()) - float(r_prev.mean())) / 10.0
        px_unc = max(5.0, half_thick + 2.0 * y_res + 2.0 * slope_px + 2.0)
        y = (row - by) / my
        dy = abs(px_unc / my)
        out.append({
            "x": float(x),
            "y": float(y),
            "lo": float(y - dy),
            "hi": float(y + dy),
            "px": [float(col), float(row)],
            "px_lo": float(row - px_unc),
            "px_hi": float(row + px_unc),
            "note": f"x_res={x_res:.2f}px, y_res={y_res:.2f}px, half_thick={half_thick:.1f}px",
        })
    return out
