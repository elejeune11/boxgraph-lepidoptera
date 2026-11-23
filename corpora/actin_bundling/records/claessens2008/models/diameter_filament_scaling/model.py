import math


def run(inputs, settings):
    return {"D": math.sqrt(inputs["n_f"])}


def symbolic():
    import sympy as sp

    n_f = sp.Symbol("n_f", positive=True)
    return {"D": sp.sqrt(n_f)}
