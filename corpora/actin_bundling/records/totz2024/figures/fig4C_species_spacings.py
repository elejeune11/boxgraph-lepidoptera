"""Calibrated 1200 dpi reading of the 5 vertical black dashed species actin-bundle spacings in Fig. 4C."""

import numpy as np
from boxgraph import digitize as D


def _find_ticks(strip, offset, axis):
    b = strip.mean(axis=axis) > 0.6
    out, i = [], 0
    while i < len(b):
        if b[i]:
            j = i
            while j + 1 < len(b) and b[j + 1]:
                j += 1
            if 4 <= j - i + 1 <= 25:
                out.append((i + j) / 2.0 + 0.5 + offset)
            i = j + 1
        else:
            i += 1
    return out


def read(image, at):
    gray = image.astype(float).mean(axis=2)
    dark = gray < 80
    H, W = dark.shape
    sub = dark[:, :int(W * 0.87)]
    hr = np.where(sub.mean(axis=1) > 0.6)[0]
    r_top = (int(hr[hr < H // 2][0]), int(hr[hr < H // 2][-1]))
    r_bot = (int(hr[hr > H // 2][0]), int(hr[hr > H // 2][-1]))
    vc = np.where(sub[r_top[0]:r_bot[1] + 1, :].mean(axis=0) > 0.6)[0]
    mid_c = (int(vc[0]) + int(vc[-1])) // 2
    c_left = (int(vc[vc < mid_c][0]), int(vc[vc < mid_c][-1]))
    c_right = (int(vc[vc > mid_c][0]), int(vc[vc > mid_c][-1]))

    # Calibrate x-axis from the 3 long major interior ticks (0.5, 1.0, 1.5 um) + right spine (2.0 um);
    # note the left spine sits at L = 0.01 um (Mathematica plot domain starts at 0.01 to avoid L=0 singularity).
    xb_long = _find_ticks(dark[r_bot[0] - 42:r_bot[0] - 36, c_left[1] + 10:c_right[0] - 10], c_left[1] + 10, axis=0)
    assert len(xb_long) == 3, f"expected 3 long x-ticks (0.5, 1.0, 1.5), got {xb_long}"
    c_r_mid = 0.5 * (c_right[0] + c_right[1]) + 0.5
    xa = D.Axis("linear", [(xb_long[0], 0.5), (xb_long[1], 1.0), (xb_long[2], 1.5), (c_r_mid, 2.0)])

    interior = gray[r_top[1] + 25:r_bot[0] - 25, c_left[1] + 15:c_right[0] - 15]
    col_dark = (interior < 60).mean(axis=0) > 0.25
    cols = np.nonzero(col_dark)[0] + c_left[1] + 15
    lines = []
    for c in cols:
        if not lines or c - lines[-1][-1] > 8:
            lines.append([int(c)])
        else:
            lines[-1].append(int(c))
    assert len(lines) == 5, f"expected 5 vertical dashed lines, found {len(lines)}"

    labels = [
        "L_2,A (Colias eurytheme lower flat surface)",
        "L_2,B (Agraulis vanillae flat inter-bundle spacing)",
        "L_1,A (Colias eurytheme upper ridged surface)",
        "L_1,B (Agraulis vanillae ridged inter-bundle spacing)",
        "L (Vanessa cardui reference inter-bundle spacing)",
    ]
    r_mid = 0.5 * (r_top[1] + r_bot[0])
    err_px = 3.0 + xa.residual_px
    out = []
    for idx, (cl, lbl) in enumerate(zip(lines, labels), start=1):
        cm = 0.5 * (cl[0] + cl[-1]) + 0.5
        L_val = xa.to_value(cm)
        lo, hi = sorted((xa.to_value(cm - err_px), xa.to_value(cm + err_px)))
        out.append({
            "x": float(idx),
            "y": L_val,
            "lo": lo,
            "hi": hi,
            "px": [cm, r_mid],
            "note": f"{lbl} (width {cl[-1]-cl[0]+1} px; res_x={xa.residual_px:.2f} px)",
        })
    return out
