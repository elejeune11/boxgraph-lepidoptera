import math


def run(inputs, settings):
    return {"R_tube": math.sqrt(inputs["kappa"] / (2.0 * inputs["sigma"]))}


def symbolic():
    import sympy as sp

    kappa, sigma = sp.symbols("kappa sigma", positive=True)
    return {"R_tube": sp.sqrt(kappa / (2 * sigma))}
