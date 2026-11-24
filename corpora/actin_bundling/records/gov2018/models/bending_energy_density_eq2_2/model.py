def run(inputs, settings):
    kappa = inputs["kappa"]
    kappa_p = inputs["kappa_p"]
    n = inputs["n"]
    kappa_eff = kappa * (1.0 - n) + kappa_p * n
    f_bend = 0.5 * kappa_eff * (inputs["H"] - n * inputs["H0"]) ** 2
    return {"f_bend": f_bend, "kappa_eff": kappa_eff}


def symbolic():
    import sympy as sp

    kappa, kappa_p, H, n, H0 = sp.symbols("kappa kappa_p H n H0", real=True)
    kappa_eff = kappa * (1 - n) + kappa_p * n
    return {
        "f_bend": sp.Rational(1, 2) * kappa_eff * (H - n * H0) ** 2,
        "kappa_eff": kappa_eff,
    }
