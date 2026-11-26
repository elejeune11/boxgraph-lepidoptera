import numpy as np


def _group_runs(indices, gap=5):
    groups = []
    for idx in indices:
        if not groups or idx > groups[-1][-1] + gap:
            groups.append([int(idx)])
        else:
            groups[-1].append(int(idx))
    return groups


def _calibrate_panel_bc(image):
    gray = image.astype(float).mean(axis=2)
    h, w = gray.shape
    vfrac = (gray[int(0.2 * h):int(0.75 * h), :] < 150).mean(axis=0)
    vgroups = _group_runs(np.where(vfrac > 0.80)[0], gap=6)
    x0 = float(np.mean(vgroups[0]))
    x1 = float(np.mean(vgroups[-1]))
    hfrac = (gray[:, int(x0) + 40:int(x1) - 40] < 150).mean(axis=1)
    hgroups = _group_runs(np.where(hfrac > 0.80)[0], gap=6)
    y0 = float(np.mean(hgroups[0]))
    y1 = float(np.mean(hgroups[-1]))
    y1_top = int(hgroups[-1][0])

    # x ticks: 0, 5, 10, 15 (frame ends + bottom inward ticks)
    xtick_strip = gray[int(y1) - 28:int(y1) - 10, int(x0) + 100:int(x1) - 100]
    xdark = (xtick_strip < 150).mean(axis=0)
    x_pts = [(x0, 0.0), (x1, 15.0)]
    for val in (5.0, 10.0):
        pred = x0 + (val / 15.0) * (x1 - x0)
        lo_i = max(0, int(pred - 60) - (int(x0) + 100))
        hi_i = min(xtick_strip.shape[1], int(pred + 60) - (int(x0) + 100))
        cols = np.where(xdark[lo_i:hi_i] > 0.55)[0]
        if len(cols):
            runs = _group_runs(cols, gap=4)
            best = min(runs, key=lambda g: abs((np.mean(g) + lo_i + int(x0) + 100) - pred))
            x_pts.append((float(np.mean(best) + lo_i + int(x0) + 100), val))
    xp = np.array([p[0] for p in x_pts])
    xv = np.array([p[1] for p in x_pts])
    mx, bx = np.polyfit(xv, xp, 1)
    x_res_px = float(np.max(np.abs(xp - (mx * xv + bx))))

    # y ticks: 0 (y1), 10, 20, 30, 40 (left inward ticks; top border is 45)
    ytick_strip = gray[int(y0) + 40:int(y1) - 40, int(x0) + 10:int(x0) + 28]
    ydark = (ytick_strip < 150).mean(axis=1)
    y_pts = [(y1, 0.0)]
    for val in (10.0, 20.0, 30.0, 40.0):
        pred = y1 - (val / 45.0) * (y1 - y0)
        lo_i = max(0, int(pred - 50) - (int(y0) + 40))
        hi_i = min(ytick_strip.shape[0], int(pred + 50) - (int(y0) + 40))
        rows = np.where(ydark[lo_i:hi_i] > 0.55)[0]
        if len(rows):
            runs = _group_runs(rows, gap=4)
            best = min(runs, key=lambda g: abs((np.mean(g) + lo_i + int(y0) + 40) - pred))
            y_pts.append((float(np.mean(best) + lo_i + int(y0) + 40), val))
    yp = np.array([p[0] for p in y_pts])
    yv = np.array([p[1] for p in y_pts])
    my, by = np.polyfit(yv, yp, 1)
    y_res_px = float(np.max(np.abs(yp - (my * yv + by))))
    return gray, (x0, x1, y0, y1, y1_top), (mx, bx, x_res_px), (my, by, y_res_px)


def read(image, at):
    gray, (x0, x1, y0, y1, y1_top), (mx, bx, x_res), (my, by, y_res) = _calibrate_panel_bc(image)
    out = []
    for xval in at:
        col = int(round(mx * xval + bx))
        # Search for the solid I/(I+B) curve in cl in [0.15, 2.5], staying clearly above bottom axis
        r_hi = min(y1_top, int(round(my * 0.12 + by)))
        r_lo = int(round(my * 2.5 + by))
        patch = gray[r_lo:r_hi, col - 2:col + 3].mean(axis=1)
        dark_rows = np.where(patch < 155)[0] + r_lo
        if len(dark_rows) == 0:
            out.append({"x": float(xval), "y": None, "note": "I/(I+B) curve not found"})
            continue
        runs = _group_runs(dark_rows, gap=4)
        # Pick the lowest dark curve (largest row index) corresponding to the I/(I+B) boundary,
        # avoiding the arrow above it at cr=2.5-3.0
        y_approx = 0.13 * (xval ** 2)
        row_pred = my * y_approx + by
        best = min(runs, key=lambda g: abs(np.mean(g) - row_pred))
        row_c = float(np.mean(best))
        half_w = 0.5 * (best[-1] - best[0] + 1)
        err_px = max(2.5, half_w + y_res + 1.0)
        yval = (row_c - by) / my
        y_a = (row_c - err_px - by) / my
        y_b = (row_c + err_px - by) / my
        out.append({
            "x": float(xval),
            "y": float(yval),
            "lo": float(min(y_a, y_b)),
            "hi": float(max(y_a, y_b)),
            "px": [float(col), row_c],
            "px_lo": float(row_c - err_px),
            "px_hi": float(row_c + err_px),
            "note": f"line_thick={best[-1]-best[0]+1}px, y_cal_res={y_res:.1f}px",
        })
    return out
