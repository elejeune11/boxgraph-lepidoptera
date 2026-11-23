"""Calibrated reader for Fig. 2(b), red dotted-circle curve (Gamma_A = 4)."""
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

    # Top x-ticks (0.5, 1.5) + left (0.0) and right (2.0) borders
    top_t = _groups(np.where((gray[int(y0) + 6:int(y0) + 18, int(x0) + 25:int(x1) - 25] < 100).mean(axis=0) > 0.5)[0] + int(x0) + 25)
    if len(top_t) != 2:
        raise RuntimeError(f"Expected 2 interior x-ticks in Fig. 2(b), found {len(top_t)}")
    x_px = np.array([x0, float(top_t[0].mean()), float(top_t[1].mean()), x1])
    x_val = np.array([0.0, 0.5, 1.5, 2.0])
    mx, bx = np.polyfit(x_val, x_px, 1)
    x_res = float(np.max(np.abs(x_px - (mx * x_val + bx))))

    # Left y-ticks: 6, 4, 2, plus F=1 line, plus bottom border 0
    lft_t = _groups(np.where((gray[int(y0) + 25:int(y1) - 25, int(x0) + 6:int(x0) + 18] < 100).mean(axis=1) > 0.5)[0] + int(y0) + 25)
    if len(lft_t) != 4:
        raise RuntimeError(f"Expected 4 left y-marks (6,4,2,F=1) in Fig. 2(b), found {len(lft_t)}")
    y_px = np.array([float(g.mean()) for g in lft_t] + [y1])
    y_val = np.array([6.0, 4.0, 2.0, 1.0, 0.0])
    my, by = np.polyfit(y_val, y_px, 1)
    y_res = float(np.max(np.abs(y_px - (my * y_val + by))))
    return (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res)


def read(image, at):
    (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res) = _calibrate(image)
    red = (image[:, :, 0] > 180) & (image[:, :, 1] < 80) & (image[:, :, 2] < 80)
    mask = np.zeros_like(red)
    mask[int(y0) + 4:int(y1) - 4, int(x0) + 4:int(x1) - 4] = True
    mask[:int(y0 + 0.36 * (y1 - y0)), int(x0 + 0.55 * (x1 - x0)):] = False
    red = red & mask

    out = []
    for x in at:
        col = int(round(mx * x + bx))
        win = red[:, max(int(x0) + 4, col - 8):min(int(x1) - 4, col + 9)]
        rr, _ = np.nonzero(win)
        if len(rr) == 0:
            win = red[:, max(int(x0) + 4, col - 35):min(int(x1) - 4, col + 36)]
            rr, _ = np.nonzero(win)
            if len(rr) == 0:
                out.append({"x": float(x), "y": None, "note": "no red pixels found"})
                continue
        row = float(np.median(rr))
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
