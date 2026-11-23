"""Anisotropic circumferential friction gamma_y^f = (gamma_y^f/gamma_x^f) * gamma_x^f (Eqs. S1a-S1b)."""

import math


def run(inputs, settings):
    gamma_x = 3.0 * math.pi * inputs["eta_cyto"] * (3.0 * inputs["d_f"] + 0.4 * inputs["L0_f"])
    return {"gamma_y_f": inputs["friction_asymmetry"] * gamma_x}


def symbolic():
    import sympy as sp

    eta_cyto, d_f, L0_f, friction_asymmetry = sp.symbols(
        "eta_cyto d_f L0_f friction_asymmetry", positive=True
    )
    gamma_x = 3 * sp.pi * eta_cyto * (3 * d_f + sp.Rational(2, 5) * L0_f)
    return {"gamma_y_f": friction_asymmetry * gamma_x}
