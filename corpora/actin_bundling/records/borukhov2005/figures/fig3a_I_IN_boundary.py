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
    vfrac = (gray[int(0.2 * h):int(0.75 * h), :] < 150).mean(axis=0)
    vgroups = _group_runs(np.where(vfrac > 0.80)[0], gap=6)
    x0 = float(np.mean(vgroups[0]))
    x1 = float(np.mean(vgroups[-1]))
    hfrac = (gray[:, int(x0) + 40:int(x1) - 40] < 150).mean(axis=1)
    hgroups = _group_runs(np.where(hfrac > 0.80)[0], gap=6)
    y0 = float(np.mean(hgroups[0]))
    y1 = float(np.mean(hgroups[-1]))

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
    # Index x is c_linkers (vertical axis on page); value y is c_rods (horizontal axis on page).
    gray, (x0, x1, y0, y1), (mx, bx, x_res), (my, by, y_res) = _calibrate_panel_a(image)
    out = []
    for cl_val in at:
        row = int(round(my * cl_val + by))
        # Search around c_rods in [3.5, 5.0] for the solid vertical I/(I+N) line
        c_lo = int(mx * 3.5 + bx)
        c_hi = int(mx * 5.0 + bx)
        patch = gray[max(int(y0) + 5, row - 3):min(int(y1) - 5, row + 4), c_lo:c_hi].mean(axis=0)
        dark_cols = np.where(patch < 150)[0] + c_lo
        if len(dark_cols) == 0:
            out.append({"x": float(cl_val), "y": None, "note": "I/(I+N) vertical line not found"})
            continue
        runs = _group_runs(dark_cols, gap=4)
        pred_col = mx * 4.27 + bx
        best = min(runs, key=lambda g: abs(np.mean(g) - pred_col))
        col_c = float(np.mean(best))
        half_w = 0.5 * (best[-1] - best[0] + 1)
        err_px = max(2.0, half_w + x_res + 1.0)
        cr_val = (col_c - bx) / mx
        cr_lo = (col_c - err_px - bx) / mx
        cr_hi = (col_c + err_px - bx) / mx
        out.append({
            "x": float(cl_val),
            "y": float(cr_val),
            "lo": float(min(cr_lo, cr_hi)),
            "hi": float(max(cr_lo, cr_hi)),
            "px": [col_c, float(row)],
            "note": f"line_thick={best[-1]-best[0]+1}px, x_cal_res={x_res:.1f}px",
        })
    return out
