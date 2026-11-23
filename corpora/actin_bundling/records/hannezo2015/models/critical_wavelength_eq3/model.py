import math


def run(inputs, settings):
    D_si = inputs["D"] * 1e-12
    lam_m = 2.0 * math.pi * ((D_si * inputs["tau"] * inputs["eta"] / inputs["xi"]) ** 0.25)
    return {"lambda_c": lam_m * 1e6}


def symbolic():
    import sympy as sp

    D, tau, eta, xi = sp.symbols("D tau eta xi", positive=True)
    D_si = D * sp.Rational(1, 10**12)
    return {"lambda_c": 2 * sp.pi * (D_si * tau * eta / xi) ** sp.Rational(1, 4) * 10**6}
