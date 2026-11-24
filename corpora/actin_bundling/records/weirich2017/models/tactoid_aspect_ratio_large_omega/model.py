import math


def run(inputs, settings):
    return {"aspect_ratio": 2.0 * math.sqrt(inputs["omega"])}


def symbolic():
    import sympy as sp

    omega = sp.Symbol("omega", positive=True)
    return {"aspect_ratio": 2 * sp.sqrt(omega)}
