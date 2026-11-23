import math


def run(inputs, settings):
    return {"a_H": 4.0 * math.pi / (math.sqrt(3.0) * inputs["q_10"])}


def symbolic():
    import sympy as sp

    q_10 = sp.Symbol("q_10", positive=True)
    return {"a_H": 4 * sp.pi / (sp.sqrt(3) * q_10)}
