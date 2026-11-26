import math


def run(inputs, settings):
    c2 = math.cos(inputs["theta"]) ** 2
    s2 = 1.0 - c2
    q_xx = 1.5 * (c2 - 1.0 / 3.0)
    q_perp = 1.5 * (0.5 * s2 - 1.0 / 3.0)
    return {"S": max(q_xx, q_perp)}


def symbolic():
    import sympy as sp
    N, theta = sp.symbols("N theta", positive=True)
    return {"S": sp.Rational(3, 2) * (sp.cos(theta) ** 2 - sp.Rational(1, 3)) + 0 * N}
