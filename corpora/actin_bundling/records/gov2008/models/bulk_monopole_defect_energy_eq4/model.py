import math


def run(inputs, settings):
    R = inputs["R"]
    a = inputs["a"]
    k_b_prime = inputs["k_b_prime"]
    rho = inputs["rho"]
    x = R / a
    return {"E_defects_b": 2.0 * math.pi * k_b_prime * (1.0 / (x * x)) * (R ** 4) * (rho ** 2) * math.log(x)}


def symbolic():
    import sympy as sp

    R, a, k_b_prime, rho = sp.symbols("R a k_b_prime rho", positive=True)
    return {"E_defects_b": 2 * sp.pi * k_b_prime * (1 / (R / a) ** 2) * R**4 * rho**2 * sp.log(R / a)}
