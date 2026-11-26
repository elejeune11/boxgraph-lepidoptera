def run(inputs, settings):
    gx, gy = inputs["t_gillespie"]["x"], inputs["t_gillespie"]["y"]
    nx, ny = inputs["t_nrm"]["x"], inputs["t_nrm"]["y"]
    nmap = {round(x, 6): y for x, y in zip(nx, ny)}
    ratios = []
    for x, y in zip(gx, gy):
        if x > 0.0 and round(x, 6) in nmap and nmap[round(x, 6)] > 0.0:
            ratios.append((x, y / nmap[round(x, 6)]))
    ratios.sort()
    return {
        "speedup_t10": float(ratios[-1][1]),
        "speedup_max": float(max(r for _, r in ratios)),
    }
