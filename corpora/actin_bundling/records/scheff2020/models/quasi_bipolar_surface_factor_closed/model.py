import math


def run(inputs, settings):
    x = float(inputs["x"])
    if x > 1e3:
        inv = 1.0 / x
        inv2 = inv * inv
        f_s = math.pi * inv * (8.0 / 3.0 + (8.0 / 5.0) * inv2 - (8.0 / 21.0) * inv2 * inv2)
        return {"f_s": f_s}
    x2 = x * x
    f_s = (2.0 * math.pi / x2) * (x + x * x2 - (x2 * x2 - 1.0) * math.atan(1.0 / x))
    return {"f_s": f_s}


def symbolic():
    import sympy as sp

    x = sp.Symbol("x", positive=True)
    return {"f_s": (2 * sp.pi / x**2) * (x + x**3 - (x**4 - 1) * sp.atan(1 / x))}
