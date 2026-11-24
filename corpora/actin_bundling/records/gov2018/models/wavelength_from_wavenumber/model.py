import math


def run(inputs, settings):
    return {"lam": 2.0 * math.pi / inputs["q"]}


def symbolic():
    import sympy as sp

    q = sp.Symbol("q", positive=True)
    return {"lam": 2 * sp.pi / q}
