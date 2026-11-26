import math


def run(inputs, settings):
    R = inputs["R"]
    a = inputs["a"]
    eps_b = inputs["eps_b"]
    rho = inputs["rho"]
    lam = inputs["lam"]
    E_c = inputs["E_c"]
    beta = inputs["beta"]
    D_R = inputs["D_R"]
    E_0 = inputs["E_0"]
    x = R / a
    e_total = -math.pi * x * eps_b * (x - 2.0) / (x * x) + rho * math.pi * ((R - lam) ** 2) * E_c / (x * x)
    return {"P": math.exp(-2.0 * beta * (e_total - E_0) / D_R)}


def symbolic():
    import sympy as sp

    R, a, eps_b, rho, lam, E_c, beta, D_R, E_0 = sp.symbols("R a eps_b rho lam E_c beta D_R E_0", positive=True)
    x = R / a
    e_total = -sp.pi * x * eps_b * (x - 2) / x**2 + rho * sp.pi * (R - lam) ** 2 * E_c / x**2
    return {"P": sp.exp(-2 * beta * (e_total - E_0) / D_R)}
