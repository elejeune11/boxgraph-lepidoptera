import math


def run(inputs, settings):
    phi = float(inputs["phi"])
    beta = float(inputs["beta"])
    eps_kBT = float(inputs["eps_kBT"])
    w = beta * phi * math.exp(eps_kBT)
    return {"p_fa": w / (1.0 + w)}


def symbolic():
    import sympy as sp

    phi, beta, eps_kBT = sp.symbols("phi beta eps_kBT", positive=True)
    w = beta * phi * sp.exp(eps_kBT)
    return {"p_fa": w / (1 + w)}
