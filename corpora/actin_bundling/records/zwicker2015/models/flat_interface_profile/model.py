import math


def run(inputs, settings):
    return {"phi_B": 0.5 * inputs["psi"] * (1.0 + math.tanh(inputs["x"] / inputs["w"]))}


def symbolic():
    import sympy as sp

    psi, x, w = sp.symbols("psi x w", positive=True)
    return {"phi_B": sp.Rational(1, 2) * psi * (1 + sp.tanh(x / w))}
