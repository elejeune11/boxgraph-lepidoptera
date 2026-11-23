"""Calibrated reader for Fig. 3(a), solid red curve (Bundle size vs. Gamma_A)."""
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
    if len(top_t) != 4:
        raise RuntimeError(f"Expected 4 x-ticks (15,20,25,30) in Fig. 3(a), found {len(top_t)}")
    x_px = np.array([float(g.mean()) for g in top_t])
    x_val = np.array([15.0, 20.0, 25.0, 30.0])
    mx, bx = np.polyfit(x_val, x_px, 1)
    x_res = float(np.max(np.abs(x_px - (mx * x_val + bx))))

    rgt_t = _groups(np.where((gray[int(y0) + 25:int(y1) - 25, int(x1) - 18:int(x1) - 6] < 100).mean(axis=1) > 0.5)[0] + int(y0) + 25)
    if len(rgt_t) != 2:
        raise RuntimeError(f"Expected 2 interior y-ticks (20,10) in Fig. 3(a), found {len(rgt_t)}")
    y_px = np.array([float(rgt_t[0].mean()), float(rgt_t[1].mean()), y1])
    y_val = np.array([20.0, 10.0, 0.0])
    my, by = np.polyfit(y_val, y_px, 1)
    y_res = float(np.max(np.abs(y_px - (my * y_val + by))))
    return (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res)


def read(image, at):
    (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res) = _calibrate(image)
    red = (image[:, :, 0] > 180) & (image[:, :, 1] < 80) & (image[:, :, 2] < 80)
    mask = np.zeros_like(red)
    mask[int(y0) + 4:int(y1) - 4, int(x0) + 4:int(x1) - 4] = True
    red = red & mask

    # Horizontal-scan centerline for steep segments (avoids miter-join and end-cap vertical bias)
    rr_all, cc_all = np.nonzero(red)
    row_top = int(rr_all.min())
    col_right = int(cc_all.max())
    h_rows, h_cols = [], []
    for r in range(row_top + 12, int(y1) - 6):
        cs = np.nonzero(red[r, :])[0]
        if len(cs) > 0 and (cs.max() - cs.min()) <= 24:
            h_rows.append(float(r))
            h_cols.append(float(cs.mean()))
    h_rows = np.array(h_rows)
    h_cols = np.array(h_cols)

    out = []
    for x in at:
        col_f = mx * x + bx
        col = int(round(col_f))
        rows = np.nonzero(red[:, col])[0]
        if len(rows) == 0:
            win = red[:, max(int(x0) + 4, col - 3):min(int(x1) - 4, col + 4)]
            rows, _ = np.nonzero(win)
        if len(rows) == 0:
            out.append({"x": float(x), "y": None, "note": "no red curve pixels found"})
            continue

        span = int(rows.max() - rows.min())
        near_right_cap = (col_f >= col_right - 10.0)
        if near_right_cap and len(h_rows) >= 20:
            # Linear centerline fit on clean interior of the final steep segment (away from cap and kink)
            sel = (h_rows >= row_top + 15) & (h_rows <= row_top + 300)
            p_end = np.polyfit(h_cols[sel], h_rows[sel], 1)
            row = float(np.polyval(p_end, col_f))
            slope_px = abs(float(p_end[0]))
            half_thick = 4.5
            px_unc = max(5.0, half_thick + 2.0 * y_res + 2.0 * slope_px + 2.0)
        elif span > 16 and len(h_cols) >= 5 and h_cols.min() <= col_f <= h_cols.max():
            order = np.argsort(h_cols)
            hc, hr = h_cols[order], h_rows[order]
            row = float(np.interp(col_f, hc, hr))
            sel = np.abs(hc - col_f) <= 15.0
            slope_px = abs(float(np.polyfit(hc[sel] - col_f, hr[sel], 1)[0])) if sel.sum() >= 4 else 0.0
            half_thick = 0.5 * float(span)
            px_unc = max(5.0, 2.0 * y_res + 2.0 * slope_px + 4.5)
        else:
            row = float(rows.mean())
            half_thick = 0.5 * float(span)
            slope_px = 0.0
            r_prev = np.nonzero(red[:, max(int(x0) + 4, col - 5)])[0]
            r_next = np.nonzero(red[:, min(int(x1) - 4, col + 5)])[0]
            if len(r_prev) > 0 and len(r_next) > 0:
                slope_px = abs(float(r_next.mean()) - float(r_prev.mean())) / 10.0
            px_unc = max(4.5, half_thick + 2.0 * y_res + 2.0 * slope_px + 2.0)

        y = (row - by) / my
        dy = abs(px_unc / my)
        out.append({
            "x": float(x),
            "y": float(y),
            "lo": float(y - dy),
            "hi": float(y + dy),
            "px": [float(col_f), float(row)],
            "px_lo": float(row - px_unc),
            "px_hi": float(row + px_unc),
            "note": f"x_res={x_res:.2f}px, y_res={y_res:.2f}px, half_thick={half_thick:.1f}px",
        })
    return out
