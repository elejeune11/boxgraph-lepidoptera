import math


def run(inputs, settings):
    return {"lam": 2.0 * math.pi * inputs["R_0"] / inputs["k"]}


def symbolic():
    import sympy as sp

    R_0, k = sp.symbols("R_0 k", positive=True)
    return {"lam": 2 * sp.pi * R_0 / k}
