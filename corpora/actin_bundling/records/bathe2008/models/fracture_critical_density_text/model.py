"""Critical fiber fracture density scaling (p. 5 col 2 & SM p. 11): phi* ~ M = sqrt(N)/2."""

import math


def run(inputs, settings):
    return {"phi_star": 0.5 * math.sqrt(inputs["N"])}


def symbolic():
    import sympy as sp

    N = sp.symbols("N", positive=True)
    return {"phi_star": sp.sqrt(N) / 2}
