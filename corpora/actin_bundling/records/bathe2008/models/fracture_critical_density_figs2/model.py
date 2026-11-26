"""Critical fiber fracture density from Fig. S2 caption (p. 11): phi* = M / 8 = sqrt(N) / 16."""

import math


def run(inputs, settings):
    M = 0.5 * math.sqrt(inputs["N"])
    return {"phi_star": M / 8.0}


def symbolic():
    import sympy as sp

    N = sp.symbols("N", positive=True)
    return {"phi_star": sp.sqrt(N) / 16}
