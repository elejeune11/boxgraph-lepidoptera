import math


def run(inputs, settings):
    a = inputs["alpha"]
    phi = inputs["phi"]
    k = inputs["i"] * math.pi
    chi_i0 = 1.0 + (k * k) / (a * a) + phi * (1.0 + (a * a) / (k * k))
    return {"chi_i0": chi_i0}


def symbolic():
    import sympy as sp

    alpha, phi, i = sp.symbols("alpha phi i", positive=True)
    k = i * sp.pi
    return {"chi_i0": 1 + k**2 / alpha**2 + phi * (1 + alpha**2 / k**2)}
