import math


def run(inputs, settings):
    return {"xi": 0.3 * inputs["c"] ** (-0.5)}


def symbolic():
    import sympy as sp

    c = sp.Symbol("c", positive=True)
    return {"xi": sp.Rational(3, 10) * c ** (-sp.Rational(1, 2))}
