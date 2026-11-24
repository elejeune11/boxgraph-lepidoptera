import math


def run(inputs, settings):
    return {"w": 2.0 * math.sqrt(inputs["kappa"] / inputs["b"])}


def symbolic():
    import sympy as sp

    kappa, b = sp.symbols("kappa b", positive=True)
    return {"w": 2 * sp.sqrt(kappa / b)}
