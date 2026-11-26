import math


def run(inputs, settings):
    return {"theta": -inputs["t"] * 2.0 * math.pi / inputs["M"]}


def symbolic():
    import sympy as sp

    M, t = sp.symbols("M t", positive=True)
    return {"theta": -t * 2 * sp.pi / M}
