"""Peak amplitude S(mu) = A / (sqrt(2*pi) * sigma) of the Gaussian structure factor fit (Methods p. 11)."""

import math


def run(inputs, settings):
    return {"S_peak": inputs["A"] / (math.sqrt(2.0 * math.pi) * inputs["sigma"])}


def symbolic():
    import sympy as sp

    A, sigma = sp.symbols("A sigma", positive=True)
    return {"S_peak": A / (sp.sqrt(2 * sp.pi) * sigma)}
