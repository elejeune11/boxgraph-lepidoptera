import math


def run(inputs, settings):
    return {"L_max": 2.0 * math.pi * float(inputs["R_drop"])}


def symbolic():
    import sympy as sp

    R_drop = sp.Symbol("R_drop", positive=True)
    return {"L_max": 2 * sp.pi * R_drop}
