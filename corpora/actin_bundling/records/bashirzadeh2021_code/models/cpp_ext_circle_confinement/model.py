"""Line-by-line Python port of ext_circle::compute() in src/ext.cpp:11-22."""

import math


def run(inputs, settings):
    k = float(inputs["k_c"])
    r_pos = float(inputs["r"])
    r_c = float(inputs["r_c"])

    lsq = r_pos * r_pos
    if lsq < r_c * r_c:
        return {"U_confine": 0.0, "F_rad": 0.0}
    l = math.sqrt(lsq)
    dl = l - r_c
    energy = 0.5 * k * dl * dl
    force_rad = -k * dl  # (-k * dl / l) * pos projected onto radial unit vector
    return {"U_confine": energy, "F_rad": force_rad}


def symbolic():
    import sympy as sp

    k_c, r, r_c = sp.symbols("k_c r r_c", nonnegative=True)
    dl = r - r_c
    U = sp.Piecewise((sp.Integer(0), r < r_c), (sp.Rational(1, 2) * k_c * dl ** 2, True))
    F = sp.Piecewise((sp.Integer(0), r < r_c), (-k_c * dl, True))
    return {"U_confine": U, "F_rad": F}
