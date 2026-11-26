"""Calibrated reader for Fig. 3(a), N = 10 planar segregated upper branch (lower thick solid line)."""
import numpy as np


def _groups(idx, gap=5):
    if len(idx) == 0:
        return []
    splits = np.where(np.diff(idx) > gap)[0] + 1
    return [g for g in np.split(idx, splits)]


def _calibrate(img):
    gray = img.astype(float).mean(axis=2)
    h, w = gray.shape
    dark = gray < 128
    vg = _groups(np.where(dark.sum(axis=0) > 0.5 * h)[0], 10)
    hg = _groups(np.where(dark.sum(axis=1) > 0.5 * w)[0], 10)
    x0, x1 = float(vg[0].mean()), float(vg[-1].mean())
    y0, y1 = float(hg[0].mean()), float(hg[-1].mean())

    top_t = _groups(np.where(dark[int(y0) + 10:int(y0) + 30, int(x0) + 20:int(x1) - 20].mean(axis=0) > 0.6)[0] + int(x0) + 20)
    x_px = np.array([x0] + [float(g.mean()) for g in top_t] + [x1])
    x_val = np.linspace(4.0, 9.0, 11)
    mx, bx = np.polyfit(x_val, x_px, 1)
    x_res = float(np.max(np.abs(x_px - (mx * x_val + bx))))

    rgt_t = _groups(np.where(dark[int(y0) + 20:int(y1) - 20, int(x1) - 30:int(x1) - 10].mean(axis=1) > 0.6)[0] + int(y0) + 20)
    y_px = np.array([y0] + [float(g.mean()) for g in rgt_t] + [y1])
    y_val = np.array([2.0, 0.0, -2.0, -4.0, -6.0, -8.0, -10.0, -12.0, -14.0])
    my, by = np.polyfit(y_val, y_px, 1)
    y_res = float(np.max(np.abs(y_px - (my * y_val + by))))
    return (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res), dark


def _solid_line_rows(dark, col, r_lo, r_hi, half_len=24):
    """Detect continuous tilted solid-line strokes passing through `col`, rejecting discrete markers."""
    dx = np.arange(-half_len, half_len + 1)
    best_score = np.zeros(r_hi - r_lo, dtype=float)
    for slope in np.linspace(0.0, 0.65, 27):
        dr = np.round(slope * dx).astype(int)
        rows = np.arange(r_lo, r_hi)[:, None] + dr[None, :]
        cols = col + dx[None, :]
        hit = dark[rows, cols] | dark[rows - 1, cols] | dark[rows + 1, cols]
        score = hit.mean(axis=1)
        best_score = np.maximum(best_score, score)
    rr = np.where(best_score >= 0.94)[0] + r_lo
    return _groups(rr, 8)


def read(image, at):
    (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res), dark = _calibrate(image)
    r_lo = int(round(my * (-6.0) + by))
    r_hi = int(round(my * (-10.2) + by))

    out = []
    for x in at:
        col_f = mx * x + bx
        col = int(round(col_f))
        grps = _solid_line_rows(dark, col, r_lo, r_hi, half_len=22 if x >= 7.8 else 26)
        if not grps:
            out.append({"x": float(x), "y": None, "note": "no continuous thick solid line found"})
            continue
        # N=20 planar is grps[0] (top thick solid line); N=10 planar is immediately below it (grps[1]),
        # or the lower edge of grps[0] once the two thick lines merge at x >= 6.5.
        if len(grps) >= 2 and x < 6.4:
            g = grps[1]
            span = float(g[-1] - g[0])
            rc = float(g.mean())
            px_unc = max(8.0, 0.5 * span + 2.0 * y_res + 3.0)
        else:
            g = grps[0]
            span = float(g[-1] - g[0])
            if span > 16:
                rc = float(g[-1]) - 4.5
                px_unc = max(12.0, 2.0 * y_res + 10.0)
            else:
                rc = float(g.mean())
                px_unc = max(8.0, 0.5 * span + 2.0 * y_res + 3.0)
        y = (rc - by) / my
        dy = abs(px_unc / my)
        out.append({
            "x": float(x),
            "y": float(y),
            "lo": float(y - dy),
            "hi": float(y + dy),
            "px": [float(col_f), float(rc)],
            "px_lo": float(rc - px_unc),
            "px_hi": float(rc + px_unc),
            "note": f"x_res={x_res:.2f}px, y_res={y_res:.2f}px, span={span:.1f}px",
        })
    return out
