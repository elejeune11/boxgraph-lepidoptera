"""Circular base radius of a single-domain conical volume element partitioning a spherical vesicle of radius R into N equal-area / equal-volume cones."""

import math


def run(inputs, settings):
    r = float(inputs["R"])
    n = float(inputs["N_domains"])
    # pi * r_cone^2 = 4*pi*R^2 / N  and  (1/3)*pi*r_cone^2*R = (4/3)*pi*R^3 / N  =>  r_cone = 2*R / sqrt(N)
    return {"r_cone": 2.0 * r / math.sqrt(n)}


def symbolic():
    import sympy as sp

    R, N_domains = sp.symbols("R N_domains", positive=True)
    return {"r_cone": 2 * R / sp.sqrt(N_domains)}
