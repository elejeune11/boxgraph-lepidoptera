"""Effective Young's modulus E = kappa / (pi * (d/2)^2) from Table S1."""

import math


def run(inputs, settings):
    # 1 pN / um^2 = 1e-12 N / 1e-12 m^2 = 1 Pa
    area_um2 = math.pi * (0.5 * inputs["d"]) ** 2
    return {"E": inputs["kappa"] / area_um2}


def symbolic():
    import sympy as sp

    kappa, d = sp.symbols("kappa d", positive=True)
    return {"E": kappa / (sp.pi * (d / 2) ** 2)}
