def run(inputs, settings):
    p = inputs["p"]
    s = inputs["s"]
    n_u = max(1, int(round(inputs["n_unit"])))
    phi0 = inputs["phi0"]
    step = 360.0 * s / p
    angles = [((k - 1) * step + phi0) % 360.0 for k in range(1, n_u + 1)]
    diffs = []
    for psi in (0.0, 60.0, 120.0, 180.0, 240.0, 300.0):
        best = min(abs(((ang - psi + 180.0) % 360.0) - 180.0) for ang in angles)
        diffs.append(best)
    return {
        "max_mismatch": max(diffs),
        "total_mismatch": sum(diffs),
    }
