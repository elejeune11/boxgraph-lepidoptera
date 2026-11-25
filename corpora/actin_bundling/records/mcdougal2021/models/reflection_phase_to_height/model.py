"""Double-pass reflection interferometry phase-to-height relation (reconstructed).

h = phi * lambda_0 / (4 * pi * n), in nm for lambda_0 in nm and phi in rad.
"""

import math


def run(inputs, settings):
    phi = float(inputs["phi"])
    lam = float(inputs["lambda_0"])
    n = float(inputs["n"])
    return {"h": phi * lam / (4.0 * math.pi * n)}


def symbolic():
    import sympy as sp

    phi = sp.Symbol("phi", real=True)
    lambda_0, n = sp.symbols("lambda_0 n", positive=True)
    return {"h": phi * lambda_0 / (4 * sp.pi * n)}
