import math


def run(inputs, settings):
    return {"eta0": 2.0 * math.atan(1.0 / inputs["x"])}


def symbolic():
    import sympy as sp

    x = sp.Symbol("x", positive=True)
    return {"eta0": 2 * sp.atan(1 / x)}
