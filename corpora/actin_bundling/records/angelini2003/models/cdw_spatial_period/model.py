import math


def run(inputs, settings):
    return {"d_CDW": 2.0 * math.pi / inputs["q_z"]}


def symbolic():
    import sympy as sp

    q_z = sp.Symbol("q_z", positive=True)
    return {"d_CDW": 2 * sp.pi / q_z}
