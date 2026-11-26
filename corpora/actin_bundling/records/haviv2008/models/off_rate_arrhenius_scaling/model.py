import math


def run(inputs, settings):
    phi = float(inputs["phi"])
    alpha = float(inputs["alpha"])
    eps_kBT = float(inputs["eps_kBT"])
    return {"koff_over_D": math.exp(-alpha * phi * eps_kBT)}


def symbolic():
    import sympy as sp

    phi, alpha, eps_kBT = sp.symbols("phi alpha eps_kBT", positive=True)
    return {"koff_over_D": sp.exp(-alpha * phi * eps_kBT)}
