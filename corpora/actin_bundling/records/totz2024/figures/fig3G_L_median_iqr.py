"""Calibrated 1200 dpi reading of Fig. 3G median (gray line) and Q1-Q3 IQR bar of unscaled width L."""

import numpy as np
from scipy import ndimage
from boxgraph import digitize as D


def _find_ticks(strip, offset, axis):
    b = strip.mean(axis=axis) > 0.5
    out, i = [], 0
    while i < len(b):
        if b[i]:
            j = i
            while j + 1 < len(b) and b[j + 1]:
                j += 1
            if 5 <= j - i + 1 <= 25:
                out.append((i + j) / 2.0 + 0.5 + offset)
            i = j + 1
        else:
            i += 1
    return out


def read(image, at):
    gray = image.astype(float).mean(axis=2)
    sat = image.max(axis=2).astype(float) - image.min(axis=2).astype(float)
    dark = gray < 80
    H, W = dark.shape
    hr = np.where(dark.mean(axis=1) > 0.5)[0]
    r_top = (int(hr[hr < H // 2][0]), int(hr[hr < H // 2][-1]))
    r_bot = (int(hr[hr > H // 2][0]), int(hr[hr > H // 2][-1]))
    vc = np.where(dark[r_top[0]:r_bot[1] + 1, :].mean(axis=0) > 0.5)[0]
    mid_c = (int(vc[0]) + int(vc[-1])) // 2
    c_left = (int(vc[vc < mid_c][0]), int(vc[vc < mid_c][-1]))
    c_right = (int(vc[vc > mid_c][0]), int(vc[vc > mid_c][-1]))

    xt = _find_ticks(dark[r_top[1] + 10:r_top[1] + 30, c_left[1] + 10:c_right[0] - 10], c_left[1] + 10, axis=0)
    yl_in = _find_ticks(dark[r_top[1] + 10:r_bot[0] - 10, c_left[1] + 8:c_left[1] + 26], r_top[1] + 10, axis=1)
    assert len(xt) == 5 and len(yl_in) == 1

    r_t_mid = 0.5 * (r_top[0] + r_top[1]) + 0.5
    r_b_mid = 0.5 * (r_bot[0] + r_bot[1]) + 0.5
    ya = D.Axis("linear", [(r_t_mid, 3.0), (yl_in[0], 1.5), (r_b_mid, 0.0)])

    times = [float(t) for t in at] if at else [39.6, 40.0, 40.3, 40.7, 41.0]
    r0 = r_top[1] + 15
    r1 = r_bot[0] - 15
    err_px = 2.5 + ya.residual_px
    out = []
    for t_val, xc in zip(times, xt):
        ic = int(round(xc))
        strip = image[r0:r1, ic - 30:ic + 31]
        white = (strip[:, :, 0] > 235) & (strip[:, :, 1] > 235) & (strip[:, :, 2] > 235)
        lab, n = ndimage.label(white)
        dot_r = None
        for k in range(1, n + 1):
            rr, cc = np.nonzero(lab == k)
            if 200 <= len(rr) <= 2500:
                hb = rr.max() - rr.min() + 1
                wb = cc.max() - cc.min() + 1
                if 0.7 <= hb / wb <= 1.4 and abs(cc.mean() - 30) < 10:
                    dot_r = float(rr.mean()) + 0.5 + r0
                    break
        assert dot_r is not None

        sides = np.concatenate([image[r0:r1, ic - 110:ic - 45], image[r0:r1, ic + 45:ic + 110]], axis=1)
        sat_s = sides.max(axis=2).astype(float) - sides.min(axis=2).astype(float)
        m_s = sides.mean(axis=2)
        score = ((sat_s < 40) & (m_s > 95) & (m_s < 165)).mean(axis=1)
        dr_rel = int(round(dot_r - r0))
        lo_s, hi_s = max(0, dr_rel - 120), min(len(score), dr_rel + 120)
        best_r = lo_s + int(np.argmax(score[lo_s:hi_s]))
        med_rows = np.where(score[max(0, best_r - 10):best_r + 11] > max(0.25, score[best_r] * 0.6))[0] + max(0, best_r - 10)
        med_r = float(med_rows.mean()) + 0.5 + r0

        in_g = np.concatenate([gray[r0:r1, ic - 24:ic - 15], gray[r0:r1, ic + 15:ic + 24]], axis=1)
        in_s = np.concatenate([sat[r0:r1, ic - 24:ic - 15], sat[r0:r1, ic + 15:ic + 24]], axis=1)
        out_g = np.concatenate([gray[r0:r1, ic - 55:ic - 35], gray[r0:r1, ic + 35:ic + 55]], axis=1)
        out_bg = ndimage.median_filter(np.median(out_g, axis=1), size=25)
        in_box = (out_bg - np.median(in_g, axis=1) > 15) | ((np.median(in_s, axis=1) < 55) & (np.median(in_g, axis=1) < 172))
        in_box[max(0, dr_rel - 26):min(len(in_box), dr_rel + 27)] = True
        in_box = ndimage.binary_closing(in_box, structure=np.ones(7))
        r_up = dr_rel
        while r_up > 0 and in_box[r_up - 1]:
            r_up -= 1
        r_dn = dr_rel
        while r_dn + 1 < len(in_box) and in_box[r_dn + 1]:
            r_dn += 1
        r_q3 = float(r_up) + r0
        r_q1 = float(r_dn) + r0

        y_val = ya.to_value(med_r)
        ylo, yhi = sorted((ya.to_value(med_r - err_px), ya.to_value(med_r + err_px)))
        bar_lo, bar_hi = sorted((ya.to_value(r_q1), ya.to_value(r_q3)))
        bar_lo_out = min(ya.to_value(r_q1 + err_px), bar_lo)
        bar_hi_out = max(ya.to_value(r_q3 - err_px), bar_hi)
        out.append({
            "x": t_val,
            "y": y_val,
            "lo": ylo,
            "hi": yhi,
            "bar_lo": bar_lo,
            "bar_hi": bar_hi,
            "bar_lo_outer": bar_lo_out,
            "bar_hi_outer": bar_hi_out,
            "px": [float(xc), med_r],
            "px_lo": r_q3,
            "px_hi": r_q1,
            "note": f"gray horizontal median & Q1-Q3 IQR bar (res_y={ya.residual_px:.2f} px)",
        })
    return out
