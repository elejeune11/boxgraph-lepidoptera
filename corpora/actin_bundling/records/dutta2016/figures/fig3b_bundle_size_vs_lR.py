"""Calibrated reader for Fig. 3(b) main panel, solid red curve (Bundle size vs. lambda_R)."""
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
        raise RuntimeError(f"Expected 3 x-ticks (0.5,0.7,0.9) in Fig. 3(b), found {len(bot_t)}")
    x_px = np.array([float(g.mean()) for g in bot_t])
    x_val = np.array([0.5, 0.7, 0.9])
    mx, bx = np.polyfit(x_val, x_px, 1)
    x_res = float(np.max(np.abs(x_px - (mx * x_val + bx))))

    lft_t = _groups(np.where((gray[int(y0) + 25:int(y1) - 25, int(x0) + 6:int(x0) + 18] < 100).mean(axis=1) > 0.5)[0] + int(y0) + 25)
    if len(lft_t) != 5:
        raise RuntimeError(f"Expected 5 y-ticks (12,10,8,6,4) in Fig. 3(b), found {len(lft_t)}")
    y_px = np.array([float(g.mean()) for g in lft_t])
    y_val = np.array([12.0, 10.0, 8.0, 6.0, 4.0])
    my, by = np.polyfit(y_val, y_px, 1)
    y_res = float(np.max(np.abs(y_px - (my * y_val + by))))
    return (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res)


def read(image, at):
    (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res) = _calibrate(image)
    red = (image[:, :, 0] > 180) & (image[:, :, 1] < 80) & (image[:, :, 2] < 80)
    mask = np.zeros_like(red)
    mask[int(y0) + 4:int(y1) - 4, int(x0) + 4:int(x1) - 4] = True
    # Mask out the inset box in upper right of Fig. 3(b)
    mask[:int(y0 + 0.66 * (y1 - y0)), int(x0 + 0.25 * (x1 - x0)):] = False
    red = red & mask

    # Horizontal-scan centerline for steep left branch (avoids vertical-slice bias on near-vertical curve)
    h_rows, h_cols = [], []
    for r in range(int(y0) + 6, int(y1) - 6):
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
        if span > 16 and len(h_cols) >= 5 and h_cols.min() <= col_f <= h_cols.max():
            order = np.argsort(h_cols)
            hc, hr = h_cols[order], h_rows[order]
            row = float(np.interp(col_f, hc, hr))
            sel = np.abs(hr - row) <= 35.0
            slope_px = abs(float(np.polyfit(hc[sel] - col_f, hr[sel], 1)[0])) if sel.sum() >= 5 else 0.0
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
