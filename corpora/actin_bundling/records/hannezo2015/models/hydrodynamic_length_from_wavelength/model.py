import math


def run(inputs, settings):
    D_m2_s = inputs["D"] * 1e-12
    return {"l_h": ((inputs["d"] / (2.0 * math.pi)) ** 2) / math.sqrt(D_m2_s * inputs["tau"])}


def symbolic():
    import sympy as sp

    d, D, tau = sp.symbols("d D tau", positive=True)
    return {"l_h": (d / (2 * sp.pi)) ** 2 / sp.sqrt(D * sp.Rational(1, 10**12) * tau)}
