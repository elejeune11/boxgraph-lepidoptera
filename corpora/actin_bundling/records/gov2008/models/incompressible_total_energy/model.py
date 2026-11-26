import math


def run(inputs, settings):
    R = inputs["R"]
    a = inputs["a"]
    eps_b = inputs["eps_b"]
    K = inputs["K"]
    rho = inputs["rho"]
    L_p = inputs["L_p"]
    x = R / a
    e_bind = -math.pi * x * eps_b * (x - 2.0) / (x * x)
    e_def = K * rho * ((R / L_p) ** 4) / (x * x)
    return {"E_total": e_bind + e_def}


def symbolic():
    import sympy as sp

    R, a, eps_b, K, rho, L_p = sp.symbols("R a eps_b K rho L_p", positive=True)
    x = R / a
    return {"E_total": -sp.pi * x * eps_b * (x - 2) / x**2 + K * rho * (R / L_p) ** 4 / x**2}
