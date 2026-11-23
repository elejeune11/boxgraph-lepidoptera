"""Shared reader for Fig. 5c (Stage 15 PIV mean displacement Vx/Vy vs Time, mean +/- SD, PDF p. 9, 1200 dpi)."""

from __future__ import annotations

import numpy as np


def _split_runs(idx: np.ndarray, gap: int = 5) -> list[np.ndarray]:
    if len(idx) == 0:
        return []
    cuts = np.where(np.diff(idx) > gap)[0] + 1
    return [g for g in np.split(idx, cuts) if len(g) > 0]


def read_fig5c(image: np.ndarray, component: str):
    """Read Vx (coral red) or Vy (teal) mean line vertices and symmetric T-capped SD error bars at -40..-0* min."""
    img = image.astype(float)
    H, W, _ = img.shape
    g = img.mean(axis=2)

    hr = _split_runs(np.where((g < 150.0).sum(axis=1) > W * 0.45)[0], gap=10)[-1]
    vc = _split_runs(np.where((g < 150.0).sum(axis=0) > H * 0.45)[0], gap=10)[0]

    # Outward y-ticks (30 and 20)
    yp = (255.0 - g[:, max(0, int(vc.min()) - 18) : max(0, int(vc.min()) - 2)]).mean(axis=1)
    yp[int(hr.min()) + 5 :] = 0.0
    yt = [float(((gr + 0.5) * yp[gr]).sum() / yp[gr].sum()) for gr in _split_runs(np.where(yp > 20.0)[0], gap=4)]
    if len(yt) != 2:
        raise RuntimeError(f"Expected 2 y-ticks (30, 20) in Fig. 5c, found {len(yt)}: {yt}")

    # Outward x-ticks (-40, -30, -20, -10, 0*)
    xp = (255.0 - g[int(hr.max()) + 2 : min(H, int(hr.max()) + 18), :]).mean(axis=0)
    xp[: int(vc.max()) + 5] = 0.0
    xt = [float(((gr + 0.5) * xp[gr]).sum() / xp[gr].sum()) for gr in _split_runs(np.where(xp > 20.0)[0], gap=4)]
    x_labels = [-40.0, -30.0, -20.0, -10.0, 0.0]
    if len(xt) != len(x_labels):
        raise RuntimeError(f"Expected 5 x-ticks in Fig. 5c, found {len(xt)}: {xt}")

    bx, ax = np.polyfit(xt, x_labels, 1)
    by, ay = np.polyfit(yt, [30.0, 20.0], 1)
    res_x = float(np.max(np.abs(np.array(xt) - (np.array(x_labels) - ax) / bx)))

    R, G, B = img[:, :, 0], img[:, :, 1], img[:, :, 2]
    if component == "vx":
        line_mask = (R - G > 60.0) & (R - B > 60.0) & (R > 180.0)
        cap_mask = (R - G > 25.0) & (R - B > 25.0) & (R > 160.0)
    elif component == "vy":
        line_mask = (G - R > 60.0) & (B - R > 60.0) & (G > 120.0)
        cap_mask = (G - R > 25.0) & (B - R > 25.0) & (G > 110.0)
    else:
        raise ValueError(f"Unknown component: {component}")

    r_top_lim = 200
    r_bot_lim = int(hr.min()) - 10

    def extrap_line(cc_arr, c_target):
        pts = []
        for cc in cc_arr:
            rr = np.nonzero(line_mask[r_top_lim:r_bot_lim, cc])[0] + r_top_lim
            if len(rr):
                runs = _split_runs(rr, gap=4)
                best = max(runs, key=len)
                if len(best) >= 12:
                    pts.append((float(cc), 0.5 * (float(best[0]) + float(best[-1])) + 0.5))
        if len(pts) < 5:
            raise RuntimeError(f"Too few line columns for extrapolation at c={c_target}")
        xs_p = np.array([p[0] for p in pts])
        ys_p = np.array([p[1] for p in pts])
        m_s, b_s = np.polyfit(xs_p, ys_p, 1)
        return float(m_s * c_target + b_s)

    out = []
    for i, (xv, c_f) in enumerate(zip(x_labels, xt)):
        ci = int(round(c_f))
        if i == 0:
            r_mean = extrap_line(np.arange(ci + 22, ci + 80), c_f)
        elif i == len(xt) - 1:
            r_mean = extrap_line(np.arange(ci - 80, ci - 22), c_f)
        else:
            r_l = extrap_line(np.arange(ci - 80, ci - 22), c_f)
            r_r = extrap_line(np.arange(ci + 22, ci + 80), c_f)
            r_mean = 0.5 * (r_l + r_r)

        # Detect horizontal T-caps centered on ci (width 28-38 px vs scatter dots 24 px)
        sub_c = cap_mask[r_top_lim:r_bot_lim, ci - 22 : ci + 23]
        wing = sub_c[:, 6:15].any(axis=1) & sub_c[:, 30:39].any(axis=1) & sub_c[:, 20:25].any(axis=1)
        caps_up, caps_dn = [], []
        for gr in _split_runs(np.where(wing)[0], gap=3):
            rm = float(gr.mean()) + 0.5 + r_top_lim
            if abs(rm - r_mean) < 25.0:
                continue  # skip thick line itself
            cols_on = np.where(sub_c[gr, :].any(axis=0))[0]
            w_tot = int(cols_on[-1] - cols_on[0] + 1)
            # T-caps are wider (>= 27 px) or paired symmetrically across r_mean
            if rm < r_mean:
                caps_up.append((rm, w_tot, len(gr)))
            else:
                caps_dn.append((rm, w_tot, len(gr)))

        # Prefer a symmetric pair (d_up ≈ d_dn within 18 px) where at least one cap has width >= 27 px
        best_pair = None
        for ru, wu, _ in caps_up:
            du = r_mean - ru
            for rd, wd, _ in caps_dn:
                dd = rd - r_mean
                if abs(du - dd) <= 18.0 and max(wu, wd) >= 27:
                    score = (wu + wd) - 2.0 * abs(du - dd)
                    if best_pair is None or score > best_pair[0]:
                        best_pair = (score, ru, rd, "both-caps")
        if best_pair is not None:
            _, r_hi_val, r_lo_val, cap_note = best_pair
        else:
            # Mirror the widest single T-cap (>= 27 px)
            wide_up = [c for c in caps_up if c[1] >= 27]
            wide_dn = [c for c in caps_dn if c[1] >= 27]
            cands = wide_up + wide_dn
            if not cands:
                raise RuntimeError(f"Fig. 5c {component} t={xv}: no T-cap found")
            best_c = max(cands, key=lambda c: c[1])
            d_sd = abs(best_c[0] - r_mean)
            r_hi_val = r_mean - d_sd
            r_lo_val = r_mean + d_sd
            cap_note = "mirrored-cap"

        y_val = ay + by * r_mean
        bar_hi = ay + by * r_hi_val
        bar_lo = ay + by * r_lo_val
        err_px = 2.5
        err_y = abs(by) * err_px
        bar_pad = abs(by) * 3.5

        out.append({
            "x": xv,
            "y": y_val,
            "lo": y_val - err_y,
            "hi": y_val + err_y,
            "bar_lo": bar_lo,
            "bar_hi": bar_hi,
            "bar_lo_outer": bar_lo - bar_pad,
            "bar_hi_outer": bar_hi + bar_pad,
            "px": [c_f, r_mean],
            "px_lo": r_hi_val,
            "px_hi": r_lo_val,
            "note": f"PIV {component} mean +/- SD ({cap_note}, rx={res_x:.2f}px)",
        })
    return out
