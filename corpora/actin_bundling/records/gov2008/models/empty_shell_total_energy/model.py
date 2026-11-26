import math


def run(inputs, settings):
    R = inputs["R"]
    a = inputs["a"]
    eps_b = inputs["eps_b"]
    rho = inputs["rho"]
    lam = inputs["lam"]
    E_c = inputs["E_c"]
    x = R / a
    e_bind = -math.pi * x * eps_b * (x - 2.0) / (x * x)
    e_def = rho * math.pi * ((R - lam) ** 2) * E_c / (x * x)
    return {"E_total": e_bind + e_def}


def symbolic():
    import sympy as sp

    R, a, eps_b, rho, lam, E_c = sp.symbols("R a eps_b rho lam E_c", positive=True)
    x = R / a
    return {"E_total": -sp.pi * x * eps_b * (x - 2) / x**2 + rho * sp.pi * (R - lam) ** 2 * E_c / x**2}
