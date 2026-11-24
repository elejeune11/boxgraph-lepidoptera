import math


def run(inputs, settings):
    return {"gamma": (1.0 / 6.0) * inputs["psi"]**2 * math.sqrt(inputs["b"] * inputs["kappa"])}


def symbolic():
    import sympy as sp

    psi, b, kappa = sp.symbols("psi b kappa", positive=True)
    return {"gamma": sp.Rational(1, 6) * psi**2 * sp.sqrt(b * kappa)}
