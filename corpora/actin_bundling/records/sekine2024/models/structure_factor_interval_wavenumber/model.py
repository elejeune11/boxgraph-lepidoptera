"""Structure interval = 2*pi / mu converted from um to nm (Methods p. 11)."""

import math


def run(inputs, settings):
    return {"interval": 1000.0 * (2.0 * math.pi) / inputs["mu"]}


def symbolic():
    import sympy as sp

    mu = sp.Symbol("mu", positive=True)
    return {"interval": sp.Integer(2000) * sp.pi / mu}
