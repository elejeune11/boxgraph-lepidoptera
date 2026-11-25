import math


def run(inputs, settings):
    return {"fwhm": 2.0 * inputs["delta"] * math.sqrt(2.0 * math.log(2.0))}


def symbolic():
    import sympy as sp
    delta = sp.Symbol("delta", positive=True)
    return {"fwhm": 2 * delta * sp.sqrt(2 * sp.log(2))}
