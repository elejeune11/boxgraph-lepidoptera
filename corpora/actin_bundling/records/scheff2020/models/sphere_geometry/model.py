import math


def run(inputs, settings):
    R = float(inputs["R"])
    return {"S": 4.0 * math.pi * R * R, "V": (4.0 / 3.0) * math.pi * (R**3)}


def symbolic():
    import sympy as sp

    R = sp.Symbol("R", positive=True)
    return {"S": 4 * sp.pi * R**2, "V": sp.Rational(4, 3) * sp.pi * R**3}
