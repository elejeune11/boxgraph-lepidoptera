import math

def _ols_slope(curve, t_cut):
    lx, ly = [], []
    for x, y in zip(curve["x"], curve["y"]):
        if 0.0 < x <= t_cut and y > 0.0:
            lx.append(math.log(x))
            ly.append(math.log(y))
    n = len(lx)
    mx = sum(lx) / n
    my = sum(ly) / n
    sxx = sum((a - mx) ** 2 for a in lx)
    sxy = sum((a - mx) * (b - my) for a, b in zip(lx, ly))
    return sxy / sxx

def run(inputs, settings):
    t_cut = float(inputs["t_cut"])
    pairs = [
        (0.125, "msd_0p125"),
        (0.25,  "msd_0p25"),
        (0.5,   "msd_0p5"),
        (1.0,   "msd_1"),
        (2.0,   "msd_2"),
        (4.0,   "msd_4"),
        (8.0,   "msd_8"),
    ]
    xs = [chi for chi, _ in pairs]
    ys = [_ols_slope(inputs[k], t_cut) for _, k in pairs]
    return {"nu_vs_chi": {"x": xs, "y": ys}}
