import math


def run(inputs, settings):
    return {"k": 2.0 * math.pi * inputs["R_0"] / inputs["lam"]}


def symbolic():
    import sympy as sp

    R_0, lam = sp.symbols("R_0 lam", positive=True)
    return {"k": 2 * sp.pi * R_0 / lam}
