"""Calibrated reader for Fig. 3(b), N = 20 planar segregated branch (top thick solid line)."""
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

    top_t = _groups(np.where(dark[int(y0) + 12:int(y0) + 32, int(x0) + 50:int(x1) - 50].mean(axis=0) > 0.6)[0] + int(x0) + 50)
    if len(top_t) != 1:
        raise RuntimeError(f"Expected 1 major top x-tick (0.1) in Fig. 3(b), found {len(top_t)}")
    x_px = np.array([x0, float(top_t[0].mean()), x1])
    x_log = np.array([-2.0, -1.0, 0.0])
    mx, bx = np.polyfit(x_log, x_px, 1)
    x_res = float(np.max(np.abs(x_px - (mx * x_log + bx))))

    rgt_t = _groups(np.where(dark[int(y0) + 50:int(y1) - 50, int(x1) - 32:int(x1) - 12].mean(axis=1) > 0.6)[0] + int(y0) + 50)
    if len(rgt_t) != 1:
        raise RuntimeError(f"Expected 1 major right y-tick (0.1) in Fig. 3(b), found {len(rgt_t)}")
    y_px = np.array([y0, float(rgt_t[0].mean()), y1])
    y_log = np.array([0.0, -1.0, -2.0])
    my, by = np.polyfit(y_log, y_px, 1)
    y_res = float(np.max(np.abs(y_px - (my * y_log + by))))
    return (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res), dark


def _solid_line_rows(dark, col, r_lo, r_hi, half_len=24):
    dx = np.arange(-half_len, half_len + 1)
    best_score = np.zeros(r_hi - r_lo, dtype=float)
    # Curve rises with x -> row decreases with col (negative to slightly positive slopes)
    for slope in np.linspace(-1.4, 0.2, 41):
        dr = np.round(slope * dx).astype(int)
        rows = np.arange(r_lo, r_hi)[:, None] + dr[None, :]
        cols = col + dx[None, :]
        hit = dark[rows, cols] | dark[rows - 1, cols] | dark[rows + 1, cols]
        score = hit.mean(axis=1)
        best_score = np.maximum(best_score, score)
    rr = np.where(best_score >= 0.92)[0] + r_lo
    return _groups(rr, 10)


def read(image, at):
    (x0, y0, x1, y1), (mx, bx, x_res), (my, by, y_res), dark = _calibrate(image)
    r_lo = int(round(y0)) + 45
    r_hi = int(round(my * np.log10(0.20) + by))

    out = []
    for x in at:
        col_f = mx * np.log10(x) + bx
        col = int(round(col_f))
        grps = _solid_line_rows(dark, col, r_lo, r_hi, half_len=20 if x >= 0.39 else 25)
        if not grps:
            out.append({"x": float(x), "y": None, "note": "no top solid curve found"})
            continue
        g = grps[0]
        span = float(g[-1] - g[0])
        rc = float(g.mean())
        px_unc = max(10.0, 0.5 * span + 2.0 * y_res + 4.0)
        y = 10.0 ** ((rc - by) / my)
        # Note my < 0 so rc - px_unc gives higher y, rc + px_unc gives lower y
        y_hi = 10.0 ** (((rc - px_unc) - by) / my)
        y_lo = 10.0 ** (((rc + px_unc) - by) / my)
        out.append({
            "x": float(x),
            "y": float(y),
            "lo": float(y_lo),
            "hi": float(y_hi),
            "px": [float(col_f), float(rc)],
            "px_lo": float(rc - px_unc),
            "px_hi": float(rc + px_unc),
            "note": f"x_res={x_res:.2f}px, y_res={y_res:.2f}px, span={span:.1f}px",
        })
    return out
