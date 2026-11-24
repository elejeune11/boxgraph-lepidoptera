import math


def run(inputs, settings):
    return {"R_s": (3.0 * inputs["V_s"] / (4.0 * math.pi)) ** (1.0 / 3.0)}


def symbolic():
    import sympy as sp

    V_s = sp.symbols("V_s", positive=True)
    return {"R_s": (3 * V_s / (4 * sp.pi)) ** sp.Rational(1, 3)}
