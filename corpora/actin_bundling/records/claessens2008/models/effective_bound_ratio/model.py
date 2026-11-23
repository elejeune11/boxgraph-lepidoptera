import math


def run(inputs, settings):
    R = inputs["R"]
    b = 1.0 + R + inputs["K_d"] / inputs["c_a"]
    disc = max(0.0, b * b - 4.0 * R)
    return {"R_star": 0.5 * (b - math.sqrt(disc))}


def symbolic():
    import sympy as sp

    R, K_d, c_a = sp.symbols("R K_d c_a", positive=True)
    b = 1 + R + K_d / c_a
    return {"R_star": sp.Rational(1, 2) * (b - sp.sqrt(b**2 - 4 * R))}
