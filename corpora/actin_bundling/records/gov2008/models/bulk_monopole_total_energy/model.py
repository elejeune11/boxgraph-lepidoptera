import math


def run(inputs, settings):
    R = inputs["R"]
    a = inputs["a"]
    eps_b = inputs["eps_b"]
    k_b_prime = inputs["k_b_prime"]
    rho = inputs["rho"]
    x = R / a
    e_bind = -math.pi * x * eps_b * (x - 2.0) / (x * x)
    e_def = 2.0 * math.pi * k_b_prime * (1.0 / (x * x)) * (R ** 4) * (rho ** 2) * math.log(x)
    return {"E_total": e_bind + e_def}


def symbolic():
    import sympy as sp

    R, a, eps_b, k_b_prime, rho = sp.symbols("R a eps_b k_b_prime rho", positive=True)
    x = R / a
    return {
        "E_total": -sp.pi * x * eps_b * (x - 2) / x**2
        + 2 * sp.pi * k_b_prime * (1 / x**2) * R**4 * rho**2 * sp.log(x)
    }
