import math


def run(inputs, settings):
    ratio = inputs["delta_l_xl"] / inputs["l_g"]
    return {"theta": math.asin(ratio)}


def symbolic():
    import sympy as sp

    delta_l_xl, l_g = sp.symbols("delta_l_xl l_g", positive=True)
    return {"theta": sp.asin(delta_l_xl / l_g)}
