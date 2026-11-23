import math


def run(inputs, settings):
    gamma = inputs["gamma"]
    mu_perp_bar = inputs["mu_perp_bar"]
    return {"f_0": 3.0 * math.pi * (gamma**4 / (16.0 * mu_perp_bar)) ** (1.0 / 3.0)}


def symbolic():
    import sympy as sp
    gamma, mu_perp_bar = sp.symbols("gamma mu_perp_bar", positive=True)
    return {"f_0": 3 * sp.pi * (gamma**4 / (16 * mu_perp_bar)) ** sp.Rational(1, 3)}
