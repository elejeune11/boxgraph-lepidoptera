"""Gaussian structure-factor fit in q/(2*pi) spatial-frequency coordinates (Supplementary Table 4 & Fig. 7c)."""

import math


def run(inputs, settings):
    A = inputs["A"]
    mu_f = inputs["mu_over_2pi"]
    sigma = inputs["sigma"]
    q_f = inputs["q_over_2pi"]
    S_q = (A / (math.sqrt(2.0 * math.pi) * sigma)) * math.exp(-((q_f - mu_f) ** 2) / (2.0 * sigma ** 2))
    return {"S_q": S_q}


def symbolic():
    import sympy as sp

    A, mu_over_2pi, sigma, q_over_2pi = sp.symbols(
        "A mu_over_2pi sigma q_over_2pi", positive=True
    )
    return {
        "S_q": (A / (sp.sqrt(2 * sp.pi) * sigma))
        * sp.exp(-((q_over_2pi - mu_over_2pi) ** 2) / (2 * sigma ** 2))
    }
