import math


def run(inputs, settings):
    x = float(inputs["x"])
    if x > 1e3:
        inv = 1.0 / x
        inv2 = inv * inv
        v = math.pi * inv2 * (16.0 / 15.0 + (32.0 / 105.0) * inv2)
        return {"v": v}
    x2 = x * x
    x3 = x * x2
    x5 = x3 * x2
    v = (math.pi / (2.0 * x3)) * (
        x + (2.0 / 3.0) * x3 + x5 - (x2 - 1.0) * ((x2 + 1.0) ** 2) * math.atan(1.0 / x)
    )
    return {"v": v}


def symbolic():
    import sympy as sp

    x = sp.Symbol("x", positive=True)
    return {
        "v": (sp.pi / (2 * x**3))
        * (x + sp.Rational(2, 3) * x**3 + x**5 - (x**2 - 1) * (x**2 + 1) ** 2 * sp.atan(1 / x))
    }
