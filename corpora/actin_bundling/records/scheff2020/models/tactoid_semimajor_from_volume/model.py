import math


def _v(x):
    if x > 1e3:
        inv = 1.0 / x
        inv2 = inv * inv
        return math.pi * inv2 * (16.0 / 15.0 + (32.0 / 105.0) * inv2)
    x2 = x * x
    x3 = x * x2
    x5 = x3 * x2
    return (math.pi / (2.0 * x3)) * (
        x + (2.0 / 3.0) * x3 + x5 - (x2 - 1.0) * ((x2 + 1.0) ** 2) * math.atan(1.0 / x)
    )


def run(inputs, settings):
    return {"R": (float(inputs["V"]) / _v(float(inputs["x"]))) ** (1.0 / 3.0)}


def symbolic():
    import sympy as sp

    V, x = sp.symbols("V x", positive=True)
    v = (sp.pi / (2 * x**3)) * (
        x + sp.Rational(2, 3) * x**3 + x**5 - (x**2 - 1) * (x**2 + 1) ** 2 * sp.atan(1 / x)
    )
    return {"R": (V / v) ** sp.Rational(1, 3)}
