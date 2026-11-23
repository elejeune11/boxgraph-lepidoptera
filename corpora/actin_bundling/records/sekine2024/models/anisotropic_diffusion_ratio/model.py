"""Passive anisotropic Brownian displacement ratio <|dy|>/<|dx|> = (gamma_y^f / gamma_x^f)^(-1/2)."""

import math


def run(inputs, settings):
    return {"R_aniso": 1.0 / math.sqrt(inputs["friction_asymmetry"])}


def symbolic():
    import sympy as sp

    friction_asymmetry = sp.Symbol("friction_asymmetry", positive=True)
    return {"R_aniso": 1 / sp.sqrt(friction_asymmetry)}
