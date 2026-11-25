import math


def run(inputs, settings):
    A = float(inputs["A"])
    K = float(inputs["K"])
    return {
        "sigma_eq": math.sqrt(8.0 * K * A / 9.0),
        "xi_eq": math.sqrt(K / (2.0 * A)),
    }


def symbolic():
    import sympy as sp

    A, K = sp.symbols("A K", positive=True)
    return {
        "sigma_eq": sp.sqrt(sp.Rational(8, 9) * K * A),
        "xi_eq": sp.sqrt(K / (2 * A)),
    }
