def run(inputs, settings):
    p = inputs["p"]
    s = inputs["s"]
    n_u = max(1, int(round(inputs["n_unit"])))
    phi0 = inputs["phi0"]
    psi = inputs["psi"]
    step = 360.0 * s / p
    best = 180.0
    for k in range(1, n_u + 1):
        ang = ((k - 1) * step + phi0) % 360.0
        diff = abs(((ang - psi + 180.0) % 360.0) - 180.0)
        if diff < best:
            best = diff
    return {"delta_phi": best}
