import math


def run(inputs, settings):
    return {"c": math.sqrt(inputs["A_star"])}


def symbolic():
    import sympy as sp

    A_star = sp.Symbol("A_star", positive=True)
    return {"c": sp.sqrt(A_star)}
