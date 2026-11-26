import math


def run(inputs, settings):
    return {"V_d": (4.0 / 3.0) * math.pi * (inputs["R_d"] ** 3)}


def symbolic():
    import sympy as sp

    R_d = sp.symbols("R_d", positive=True)
    return {"V_d": sp.Rational(4, 3) * sp.pi * R_d**3}
