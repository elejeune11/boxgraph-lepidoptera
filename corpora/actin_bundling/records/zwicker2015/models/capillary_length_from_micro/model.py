import math


def run(inputs, settings):
    gamma = (1.0 / 6.0) * inputs["psi"]**2 * math.sqrt(inputs["b"] * inputs["kappa"])
    beta = 2.0 / (inputs["b"] * inputs["psi"])
    return {"gamma_beta": gamma * beta}


def symbolic():
    import sympy as sp

    b, kappa, psi = sp.symbols("b kappa psi", positive=True)
    gamma = sp.Rational(1, 6) * psi**2 * sp.sqrt(b * kappa)
    beta = 2 / (b * psi)
    return {"gamma_beta": sp.simplify(gamma * beta)}
