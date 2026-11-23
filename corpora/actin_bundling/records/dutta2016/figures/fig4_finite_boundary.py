"""Calibrated reader for Fig. 4, red dashed-circle curve (finite-bundle spinodal boundary lambda_R vs. Gamma_A).
Isolates filled red circular markers via Euclidean distance transform so steep red dashes
above circles (at Gamma_A = 80, 81) do not bias marker centroids."""
import numpy as np
from scipy.ndimage import distance_transform_edt, label


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
    red = (image[:, :, 0] > 180) & (image[:, :, 1] < 80) & (image[:, :, 2] < 80)
    mask = np.zeros_like(red)
    mask[int(y0) + 4:int(y1) - 4, int(x0) + 4:int(x1) - 4] = True
    red = red & mask

    # Isolate filled circle marker interiors (radius ~16.5 px) vs. thin dashes (half-width ~5.5 px)
    cores = distance_transform_edt(red) >= 10.0
    lbl, n_circ = label(cores)
    circ_c, circ_r = [], []
    for i in range(1, n_circ + 1):
        rr, cc = np.where(lbl == i)
        if len(rr) >= 50:
            circ_c.append(float(cc.mean()))
            circ_r.append(float(rr.mean()))
    circ_c = np.array(circ_c)
    circ_r = np.array(circ_r)
    order = np.argsort(circ_c)
    circ_c, circ_r = circ_c[order], circ_r[order]

    out = []
    for x in at:
        col = mx * x + bx
        dists = np.abs(circ_c - col)
        idx = int(np.argmin(dists)) if len(dists) > 0 else -1
        if idx >= 0 and dists[idx] <= 8.0:
            row = float(circ_r[idx])
            # Local slope from adjacent circle markers
            i0 = max(0, idx - 1)
            i1 = min(len(circ_c) - 1, idx + 1)
            slope_px = abs(float(circ_r[i1] - circ_r[i0]) / float(circ_c[i1] - circ_c[i0])) if i1 > i0 else 0.0
        elif len(circ_c) >= 2 and circ_c[0] <= col <= circ_c[-1]:
            row = float(np.interp(col, circ_c, circ_r))
            near = np.argsort(dists)[:2]
            slope_px = abs(float(circ_r[near[1]] - circ_r[near[0]]) / float(circ_c[near[1]] - circ_c[near[0]]))
        else:
            out.append({"x": float(x), "y": None, "note": "no red circle marker near x"})
            continue

        px_unc = max(5.0, 2.0 * y_res + 2.5 * slope_px + 3.0)
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
            "note": f"x_res={x_res:.2f}px, y_res={y_res:.2f}px, slope={slope_px:.2f}",
        })
    return out
