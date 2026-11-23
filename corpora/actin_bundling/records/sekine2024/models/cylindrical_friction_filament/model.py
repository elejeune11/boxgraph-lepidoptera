"""Cylindrical approximation for actin filament endpoint particle friction (SI Notes p. 2)."""

import math


def run(inputs, settings):
    # 1 Pa*s * 1 um = 1 (N/m^2)*s * 1e-6 m = 1e-6 N*s/m = 1 pN*s/um
    gamma_f = 3.0 * math.pi * inputs["eta_cyto"] * (3.0 * inputs["d_f"] + 0.4 * inputs["L0_f"])
    return {"gamma_f": gamma_f}


def symbolic():
    import sympy as sp

    eta_cyto, d_f, L0_f = sp.symbols("eta_cyto d_f L0_f", positive=True)
    return {"gamma_f": 3 * sp.pi * eta_cyto * (3 * d_f + sp.Rational(2, 5) * L0_f)}
