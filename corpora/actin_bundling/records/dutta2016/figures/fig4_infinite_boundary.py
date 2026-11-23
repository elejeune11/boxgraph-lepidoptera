"""Calibrated reader for Fig. 4, blue solid-triangle curve (infinite-bundle boundary F(k=0)=1).
Samples at midpoints between downward-triangle markers (56, 60, 64, ...) on the pure solid blue line."""
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
    if len(bot_t) != 5:
        raise RuntimeError(f"Expected 5 bottom x-ticks (60,70,80,90,100) in Fig. 4, found {len(bot_t)}")
    x_px = np.array([float(g.mean()) for g in bot_t])
    x_val = np.array([60.0, 70.0, 80.0, 90.0, 100.0])
    mx, bx = np.polyfit(x_val, x_px, 1)
    x_res = float(np.max(np.abs(x_px - (mx * x_val + bx))))

    rgt_t = _groups(np.where((gray[int(y0) + 25:int(y1) - 25, int(x1) - 18:int(x1) - 6] < 100).mean(axis=1) > 0.5)[0] + int(y0) + 25)
    if len(rgt_t) != 2:
        raise RuntimeError(f"Expected 2 right y-ticks (0.4,0.2) in Fig. 4, found {len(rgt_t)}")
    y_px = np.array([float(rgt_t[0].mean()), float(rgt_t[1].mean()), y1])
    y_val = np.array([0.4, 0.2, 0.0])
    my, by = np.polyfit(y_val, y_px, 1)
    y_res = float(np.max(np.abs(y_px - (my * y_val + by))))
    return (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res)


def read(image, at):
    (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res) = _calibrate(image)
    blue = (image[:, :, 2] > 120) & (image[:, :, 0] < 100) & (image[:, :, 1] < 120)
    mask = np.zeros_like(blue)
    mask[int(y0) + 4:int(y1) - 4, int(x0) + 4:int(x1) - 4] = True
    blue = blue & mask

    out = []
    for x in at:
        col = int(round(mx * x + bx))
        rows = np.nonzero(blue[:, col])[0]
        if len(rows) == 0:
            win = blue[:, max(int(x0) + 4, col - 3):min(int(x1) - 4, col + 4)]
            rows, _ = np.nonzero(win)
        if len(rows) == 0:
            out.append({"x": float(x), "y": None, "note": "no blue boundary pixels found"})
            continue
        row = float(rows.mean())
        half_thick = 0.5 * float(rows.max() - rows.min())
        px_unc = max(4.5, half_thick + 2.0 * y_res + 2.0)
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
