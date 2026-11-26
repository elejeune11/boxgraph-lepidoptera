import math


def run(inputs, settings):
    R = inputs["R"]
    a = inputs["a"]
    k_s_prime = inputs["k_s_prime"]
    rho = inputs["rho"]
    x = R / a
    return {"E_defects_s": 2.0 * math.pi * k_s_prime * (1.0 / (x * x)) * (R * R) * rho * math.log(x)}


def symbolic():
    import sympy as sp

    R, a, k_s_prime, rho = sp.symbols("R a k_s_prime rho", positive=True)
    return {"E_defects_s": 2 * sp.pi * k_s_prime * (1 / (R / a) ** 2) * R**2 * rho * sp.log(R / a)}
