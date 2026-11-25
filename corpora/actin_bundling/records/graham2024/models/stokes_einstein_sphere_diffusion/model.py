import math


def run(inputs, settings):
    mu = 1.0 / (6.0 * math.pi * inputs["eta"] * inputs["r"])
    kBT_pN_um = inputs["kBT"] * 1e-3
    return {"D": mu * kBT_pN_um, "mu": mu}


def symbolic():
    import sympy as sp

    eta, r, kBT = sp.symbols("eta r kBT", positive=True)
    mu = 1 / (6 * sp.pi * eta * r)
    return {"D": mu * kBT * sp.Rational(1, 1000), "mu": mu}
