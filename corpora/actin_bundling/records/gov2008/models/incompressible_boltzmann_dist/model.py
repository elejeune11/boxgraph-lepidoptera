import math


def run(inputs, settings):
    R = inputs["R"]
    a = inputs["a"]
    eps_b = inputs["eps_b"]
    K = inputs["K"]
    rho = inputs["rho"]
    L_p = inputs["L_p"]
    beta = inputs["beta"]
    D_R = inputs["D_R"]
    E_0 = inputs["E_0"]
    x = R / a
    e_total = -math.pi * x * eps_b * (x - 2.0) / (x * x) + K * rho * ((R / L_p) ** 4) / (x * x)
    return {"P": math.exp(-2.0 * beta * (e_total - E_0) / D_R)}


def symbolic():
    import sympy as sp

    R, a, eps_b, K, rho, L_p, beta, D_R, E_0 = sp.symbols(
        "R a eps_b K rho L_p beta D_R E_0", positive=True
    )
    x = R / a
    e_total = -sp.pi * x * eps_b * (x - 2) / x**2 + K * rho * (R / L_p) ** 4 / x**2
    return {"P": sp.exp(-2 * beta * (e_total - E_0) / D_R)}
