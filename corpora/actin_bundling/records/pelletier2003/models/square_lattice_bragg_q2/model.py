import math


def run(inputs, settings):
    return {"q_2": inputs["q_1"] * math.sqrt(2.0)}


def symbolic():
    import sympy as sp

    q_1 = sp.Symbol("q_1", positive=True)
    return {"q_2": q_1 * sp.sqrt(2)}
