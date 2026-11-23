import math


def run(inputs, settings):
    return {"R": (3.0 + inputs["sigma"] * 3.5 / math.sqrt(inputs["n_f"])) / inputs["u"]}


def symbolic():
    import sympy as sp

    u, n_f = sp.symbols("u n_f", positive=True)
    sigma = sp.Symbol("sigma", real=True)
    return {"R": (3 + sigma * sp.Rational(7, 2) / sp.sqrt(n_f)) / u}
