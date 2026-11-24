import math


def run(inputs, settings):
    R_0 = inputs["R_0"]
    u = inputs["u"]
    lam = inputs["lam"]
    gamma_I = inputs["gamma_I"]
    w = inputs["omega"]
    k = 2.0 * math.pi * R_0 / lam
    # (um * um^2 / um) * (nN/m) = um^2 * nN/m = 1e-21 J
    dF = (math.pi * lam * (u**2) * gamma_I / (2.0 * R_0)) * ((1.0 + w) * (k**2) - 1.0) * 1e-21
    return {"delta_F": dF}


def symbolic():
    import sympy as sp

    R_0, u, lam, gamma_I, omega = sp.symbols("R_0 u lam gamma_I omega", positive=True)
    k = 2 * sp.pi * R_0 / lam
    return {
        "delta_F": (sp.pi * lam * u**2 * gamma_I / (2 * R_0))
        * ((1 + omega) * k**2 - 1)
        * sp.Rational(1, 10**21)
    }
