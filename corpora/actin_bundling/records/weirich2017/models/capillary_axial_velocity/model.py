import math


def run(inputs, settings):
    return {"v": 2.0 * math.pi * inputs["du_dt"] / inputs["k"]}


def symbolic():
    import sympy as sp

    du_dt, k = sp.symbols("du_dt k", positive=True)
    return {"v": 2 * sp.pi * du_dt / k}
