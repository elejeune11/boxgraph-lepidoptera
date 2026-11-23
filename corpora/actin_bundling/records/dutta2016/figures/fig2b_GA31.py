"""Calibrated reader for Fig. 2(b), black solid-square curve (Gamma_A = 31)."""
import numpy as np
from scipy.ndimage import binary_dilation


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
    r_ch = image[:, :, 0].astype(float)
    g_ch = image[:, :, 1].astype(float)
    b_ch = image[:, :, 2].astype(float)
    colored = ((r_ch > 150) & (b_ch < 100)) | ((b_ch > 110) & (r_ch < 90))
    # Dilate colored mask so black borders of red circles / blue triangles are excluded
    near_color = binary_dilation(colored, iterations=12)
    blk = (r_ch < 60) & (g_ch < 60) & (b_ch < 60) & (~near_color)

    # Mask outside plot box, legend, and horizontal F=1 line
    mask = np.zeros_like(blk)
    mask[int(y0) + 6:int(y1) - 6, int(x0) + 6:int(x1) - 6] = True
    mask[:int(y0 + 0.36 * (y1 - y0)), int(x0 + 0.55 * (x1 - x0)):] = False
    mask[int(y_F1) - 8:int(y_F1) + 9, :] = False
    blk = blk & mask

    out = []
    for x in at:
        col = int(round(mx * x + bx))
        strip = blk[:, max(int(x0) + 6, col - 3):min(int(x1) - 6, col + 4)]
        rows = np.where(strip.mean(axis=1) > 0.5)[0]
        if len(rows) == 0:
            out.append({"x": float(x), "y": None, "note": "no isolated black curve pixels"})
            continue
        grps = _groups(rows, gap=14)
        # Pick the thickest cluster (black square marker at 0.1 multiples is ~45 px tall)
        g = max(grps, key=len)
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
