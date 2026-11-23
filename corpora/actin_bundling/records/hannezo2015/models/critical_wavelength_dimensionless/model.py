import math


def run(inputs, settings):
    return {"lambda_hat": 2.0 * math.pi * (inputs["phi"] ** -0.25)}


def symbolic():
    import sympy as sp

    phi = sp.Symbol("phi", positive=True)
    return {"lambda_hat": 2 * sp.pi * phi ** sp.Rational(-1, 4)}
