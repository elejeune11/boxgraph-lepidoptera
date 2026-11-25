import math


def run(inputs, settings):
    s2c = sorted(zip(inputs["theta_vs_n"]["x"], inputs["theta_vs_n"]["y"]), key=lambda p: p[1])
    sS4 = sorted(zip(inputs["theta_vs_Rmax"]["x"], inputs["theta_vs_Rmax"]["y"]), key=lambda p: p[1])
    sigma = float(inputs["sigma"])  # uN/m; with R in um, sigma/R is directly in N/m^2
    dps = []
    for (n, th1), (rmax, th2) in zip(s2c, sS4):
        th = 0.5 * (th1 + th2)
        rf = rmax / n
        dp = 2.0 * sigma * (rmax - rf * math.cos(math.radians(th))) / (rmax * rmax - rf * rf)
        dps.append(dp)
    tx = list(inputs["hist_target"]["x"])
    step = (tx[-1] - tx[0]) / max(1, len(tx) - 1)
    counts = []
    for j, xj in enumerate(tx):
        lo = round(xj)  # left bin edges 20, 30, 40, 50, 60
        hi = lo + round(step)
        if j == len(tx) - 1:
            c = sum(1 for v in dps if lo <= v <= hi)
        else:
            c = sum(1 for v in dps if lo <= v < hi)
        counts.append(float(c))
    return {"freq_curve": {"x": tx, "y": counts}}
