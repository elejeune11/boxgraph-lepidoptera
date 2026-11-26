import math


def run(inputs, settings):
    R = inputs["R"]
    a = inputs["a"]
    eps_b = inputs["eps_b"]
    k_b_prime = inputs["k_b_prime"]
    rho = inputs["rho"]
    beta = inputs["beta"]
    D_R = inputs["D_R"]
    E_0 = inputs["E_0"]
    x = R / a
    e_total = -math.pi * x * eps_b * (x - 2.0) / (x * x) + 2.0 * math.pi * k_b_prime * (a * a) * (R * R) * (rho * rho) * math.log(x)
    return {"P": math.exp(-2.0 * beta * (e_total - E_0) / D_R)}


def symbolic():
    import sympy as sp

    R, a, eps_b, k_b_prime, rho, beta, D_R, E_0 = sp.symbols(
        "R a eps_b k_b_prime rho beta D_R E_0", positive=True
    )
    x = R / a
    e_total = -sp.pi * x * eps_b * (x - 2) / x**2 + 2 * sp.pi * k_b_prime * a**2 * R**2 * rho**2 * sp.log(x)
    return {"P": sp.exp(-2 * beta * (e_total - E_0) / D_R)}
