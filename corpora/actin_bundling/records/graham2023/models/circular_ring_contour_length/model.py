import math


def run(inputs, settings):
    return {"L": 2.0 * math.pi * inputs["R"]}


def symbolic():
    import sympy as sp

    R = sp.Symbol("R", positive=True)
    return {"L": 2 * sp.pi * R}
