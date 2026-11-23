"""Cylindrical cytosolic friction gamma^X = 6*pi*eta_cyto*(3*d^X + L_0^X/5) for X in {m, cl} (SI Notes p. 2)."""

import math


def run(inputs, settings):
    gamma = 6.0 * math.pi * inputs["eta_cyto"] * (3.0 * inputs["d"] + 0.2 * inputs["L0"])
    return {"gamma": gamma}


def symbolic():
    import sympy as sp

    eta_cyto, d, L0 = sp.symbols("eta_cyto d L0", positive=True)
    return {"gamma": 6 * sp.pi * eta_cyto * (3 * d + sp.Rational(1, 5) * L0)}
