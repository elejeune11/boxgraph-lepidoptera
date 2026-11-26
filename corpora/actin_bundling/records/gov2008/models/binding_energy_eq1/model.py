import math


def run(inputs, settings):
    R = inputs["R"]
    a = inputs["a"]
    eps_b = inputs["eps_b"]
    x = R / a
    return {"E_bind": -math.pi * x * eps_b * (x - 2.0) / (x * x)}


def symbolic():
    import sympy as sp

    R, a, eps_b = sp.symbols("R a eps_b", positive=True)
    x = R / a
    return {"E_bind": -sp.pi * x * eps_b * (x - 2) / x**2}
