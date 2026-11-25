def run(inputs, settings):
    ns = list(inputs["theta_vs_n"]["x"])
    tx = list(inputs["hist_target"]["x"])
    step = (tx[-1] - tx[0]) / max(1, len(tx) - 1)
    counts = []
    for j, xj in enumerate(tx):
        lo = round(xj)  # left bin edges 2, 4, 6, 8, 10
        hi = lo + round(step)
        if j == len(tx) - 1:
            c = sum(1 for v in ns if lo <= v <= hi)
        else:
            c = sum(1 for v in ns if lo <= v < hi)
        counts.append(float(c))
    return {"freq_curve": {"x": tx, "y": counts}}
