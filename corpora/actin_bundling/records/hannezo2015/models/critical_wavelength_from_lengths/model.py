import math


def run(inputs, settings):
    return {"lambda_c": 2.0 * math.pi * math.sqrt(inputs["l_t"] * inputs["l_h"])}


def symbolic():
    import sympy as sp

    l_t, l_h = sp.symbols("l_t l_h", positive=True)
    return {"lambda_c": 2 * sp.pi * sp.sqrt(l_t * l_h)}
