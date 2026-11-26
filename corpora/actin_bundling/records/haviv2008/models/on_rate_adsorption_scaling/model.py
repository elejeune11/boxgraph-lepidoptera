import math


def run(inputs, settings):
    n_a = float(inputs["n_a"])
    phi = float(inputs["phi"])
    beta = float(inputs["beta"])
    eps_kBT = float(inputs["eps_kBT"])
    w = beta * phi * math.exp(eps_kBT)
    p_fa = w / (1.0 + w)
    return {"two_kon_over_D": n_a * p_fa}


def symbolic():
    import sympy as sp

    n_a, phi, beta, eps_kBT = sp.symbols("n_a phi beta eps_kBT", positive=True)
    w = beta * phi * sp.exp(eps_kBT)
    return {"two_kon_over_D": n_a * w / (1 + w)}
