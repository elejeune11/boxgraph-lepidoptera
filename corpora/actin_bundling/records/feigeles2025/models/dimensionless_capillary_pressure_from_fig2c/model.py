import math


def run(inputs, settings):
    nx = inputs["theta_vs_n"]["x"]
    th = inputs["theta_vs_n"]["y"]
    tx = inputs["P_star_target"]["x"]
    pairs = list(zip(nx, th))
    y_out = []
    for nt in tx:
        close = [(n, t) for n, t in pairs if abs(n - nt) <= 0.06]
        if not close:
            close = [min(pairs, key=lambda p: abs(p[0] - nt))]
        ps = [2.0 * (n - math.cos(math.radians(t))) / (n * n - 1.0) for n, t in close]
        y_out.append(sum(ps) / len(ps))
    return {"P_star_curve": {"x": list(tx), "y": y_out}}
