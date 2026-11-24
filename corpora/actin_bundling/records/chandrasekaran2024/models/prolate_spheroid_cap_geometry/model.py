def run(inputs, settings):
    # f_caps in percent -> fraction in (0, 1)
    f_total = inputs["f_caps"] / 100.0
    tol = float(settings["tol"])

    # For two symmetric polar caps |z| >= u*a of a prolate spheroid:
    # V_caps / V_total = 0.5 * (1 - u)^2 * (2 + u), strictly decreasing on u in [0, 1].
    lo, hi = 0.0, 1.0
    while (hi - lo) > tol:
        mid = 0.5 * (lo + hi)
        val = 0.5 * ((1.0 - mid) ** 2) * (2.0 + mid)
        if val > f_total:
            lo = mid
        else:
            hi = mid
    u_cap = 0.5 * (lo + hi)

    return {
        "f_single": 0.5 * inputs["f_caps"],
        "u_cap": u_cap,
        "C_thresh": inputs["C_bulk"] / 3.0,
        "C_thresh_pct": 100.0 / 3.0,
    }
