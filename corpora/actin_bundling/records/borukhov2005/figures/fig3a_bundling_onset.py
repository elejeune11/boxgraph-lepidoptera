import numpy as np


def _group_runs(indices, gap=5):
    groups = []
    for idx in indices:
        if not groups or idx > groups[-1][-1] + gap:
            groups.append([int(idx)])
        else:
            groups[-1].append(int(idx))
    return groups


def _calibrate_panel_a(image):
    gray = image.astype(float).mean(axis=2)
    h, w = gray.shape
    # Detect vertical frame borders (>80% dark across middle rows)
    vfrac = (gray[int(0.2 * h):int(0.75 * h), :] < 150).mean(axis=0)
    vgroups = _group_runs(np.where(vfrac > 0.80)[0], gap=6)
    x0 = float(np.mean(vgroups[0]))
    x1 = float(np.mean(vgroups[-1]))
    # Detect horizontal frame borders (>80% dark across interior cols)
    hfrac = (gray[:, int(x0) + 40:int(x1) - 40] < 150).mean(axis=1)
    hgroups = _group_runs(np.where(hfrac > 0.80)[0], gap=6)
    y0 = float(np.mean(hgroups[0]))
    y1 = float(np.mean(hgroups[-1]))

    # Calibrate x-axis from frame ends (0, 15) + bottom inward tick at 5 and 10
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

    # Calibrate y-axis from frame ends (0.6 at y0, 0.0 at y1) + left inward ticks (0.1..0.5)
    ytick_strip = gray[int(y0) + 40:int(y1) - 40, int(x0) + 10:int(x0) + 28]
    ydark = (ytick_strip < 150).mean(axis=1)
    y_pts = [(y0, 0.6), (y1, 0.0)]
    for val in (0.1, 0.2, 0.3, 0.4, 0.5):
        pred = y1 - (val / 0.6) * (y1 - y0)
        lo_i = max(0, int(pred - 40) - (int(y0) + 40))
        hi_i = min(ytick_strip.shape[0], int(pred + 40) - (int(y0) + 40))
        rows = np.where(ydark[lo_i:hi_i] > 0.55)[0]
        if len(rows):
            runs = _group_runs(rows, gap=4)
            best = min(runs, key=lambda g: abs((np.mean(g) + lo_i + int(y0) + 40) - pred))
            y_pts.append((float(np.mean(best) + lo_i + int(y0) + 40), val))
    yp = np.array([p[0] for p in y_pts])
    yv = np.array([p[1] for p in y_pts])
    my, by = np.polyfit(yv, yp, 1)
    y_res_px = float(np.max(np.abs(yp - (my * yv + by))))

    return gray, (x0, x1, y0, y1), (mx, bx, x_res_px), (my, by, y_res_px)


def read(image, at):
    gray, (x0, x1, y0, y1), (mx, bx, x_res), (my, by, y_res) = _calibrate_panel_a(image)
    out = []
    for xval in at:
        col = int(round(mx * xval + bx))
        # Expected linear band for the diagonal bundling onset line in Fig. 3a:
        # passes roughly from (0, 0.037) to (15, 0.402)
        y_approx = 0.037 + 0.0243 * xval
        row_pred = my * y_approx + by
        r_lo = max(int(y0) + 15, int(row_pred - 90))
        r_hi = min(int(y1) - 15, int(row_pred + 90))
        # Average 3 neighboring columns to suppress pixel jitter on diagonal line
        patch = gray[r_lo:r_hi, max(int(x0) + 5, col - 1):min(int(x1) - 5, col + 2)].mean(axis=1)
        dark_rows = np.where(patch < 150)[0] + r_lo
        if len(dark_rows) == 0:
            out.append({"x": float(xval), "y": None, "note": "bundling onset line not found"})
            continue
        runs = _group_runs(dark_rows, gap=4)
        best = min(runs, key=lambda g: abs(np.mean(g) - row_pred))
        row_c = float(np.mean(best))
        half_w = 0.5 * (best[-1] - best[0] + 1)
        err_px = max(2.0, half_w + y_res + 1.0)
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
