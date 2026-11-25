import math


def _fs(x):
    if x > 1e3:
        inv = 1.0 / x
        inv2 = inv * inv
        return math.pi * inv * (8.0 / 3.0 + (8.0 / 5.0) * inv2 - (8.0 / 21.0) * inv2 * inv2)
    x2 = x * x
    return (2.0 * math.pi / x2) * (x + x * x2 - (x2 * x2 - 1.0) * math.atan(1.0 / x))


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


def _a_spindle(x):
    if x > 1e3:
        inv = 1.0 / x
        inv2 = inv * inv
        return inv * (8.0 / 3.0 + (8.0 / 15.0) * inv2 - (8.0 / 105.0) * inv2 * inv2)
    eta0 = 2.0 * math.atan(1.0 / x)
    s0 = math.sin(eta0)
    return 2.0 * (eta0 - s0 * math.cos(eta0)) / (s0 * s0)


def run(inputs, settings):
    R = float(inputs["R"])
    x = float(inputs["x"])
    return {
        "S": R * R * _fs(x),
        "V": (R**3) * _v(x),
        "r": R / x,
        "A_spindle": R * R * _a_spindle(x),
        "A_ellipse": math.pi * R * R / x,
    }


def symbolic():
    import sympy as sp

    R, x = sp.symbols("R x", positive=True)
    fs = (2 * sp.pi / x**2) * (x + x**3 - (x**4 - 1) * sp.atan(1 / x))
    v = (sp.pi / (2 * x**3)) * (
        x + sp.Rational(2, 3) * x**3 + x**5 - (x**2 - 1) * (x**2 + 1) ** 2 * sp.atan(1 / x)
    )
    eta0 = 2 * sp.atan(1 / x)
    a_sp = 2 * (eta0 - sp.sin(eta0) * sp.cos(eta0)) / sp.sin(eta0) ** 2
    return {
        "S": R**2 * fs,
        "V": R**3 * v,
        "r": R / x,
        "A_spindle": R**2 * a_sp,
        "A_ellipse": sp.pi * R**2 / x,
    }
