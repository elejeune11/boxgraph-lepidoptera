import math


def run(inputs, settings):
    return {"d": 2.0 * math.pi / inputs["q_1"]}


def symbolic():
    import sympy as sp

    q_1 = sp.Symbol("q_1", positive=True)
    return {"d": 2 * sp.pi / q_1}
