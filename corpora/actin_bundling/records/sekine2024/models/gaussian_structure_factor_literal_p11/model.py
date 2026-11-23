"""Literal Methods p. 11 Gaussian structure-factor fit in angular wavenumber q = 2*pi*(q/(2*pi))."""

import math


def run(inputs, settings):
    A = inputs["A"]
    mu = 2.0 * math.pi * inputs["mu_over_2pi"]
    sigma = inputs["sigma"]
    q = 2.0 * math.pi * inputs["q_over_2pi"]
    S_q = (A / (math.sqrt(2.0 * math.pi) * sigma)) * math.exp(-((q - mu) ** 2) / (2.0 * sigma ** 2))
    return {"S_q": S_q}


def symbolic():
    import sympy as sp

    A, mu_over_2pi, sigma, q_over_2pi = sp.symbols(
        "A mu_over_2pi sigma q_over_2pi", positive=True
    )
    q = 2 * sp.pi * q_over_2pi
    mu = 2 * sp.pi * mu_over_2pi
    return {
        "S_q": (A / (sp.sqrt(2 * sp.pi) * sigma))
        * sp.exp(-((q - mu) ** 2) / (2 * sigma ** 2))
    }
