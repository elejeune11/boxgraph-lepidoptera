import math


def _f_and_g(a):
    if a < 1e-2:
        a2 = a * a
        f = (a**3) * (1.0 / 3.0 - a2 / 30.0 + (a2 * a2) / 840.0 - (a2**3) / 45360.0)
        g = (a**5) * (2.0 / 15.0 - 11.0 * a2 / 315.0 + 17.0 * (a2 * a2) / 3780.0 - 461.0 * (a2**3) / 1247400.0)
        return f, g
    s = math.sin(a)
    f = s - a * math.cos(a)
    g = f - (1.0 / 3.0) * (s**3)
    return f, g


def run(inputs, settings):
    a = inputs["alpha"]
    w = inputs["omega"]
    f, g = _f_and_g(a)
    F_hat = (4.0 * math.pi * f + 2.0 * math.pi * w * g) / ((2.0 * math.pi * g) ** (2.0 / 3.0))
    return {"F_hat": F_hat}


def symbolic():
    import sympy as sp

    alpha, omega = sp.symbols("alpha omega", positive=True)
    f = sp.sin(alpha) - alpha * sp.cos(alpha)
    g = f - sp.Rational(1, 3) * sp.sin(alpha) ** 3
    return {
        "F_hat": (4 * sp.pi * f + 2 * sp.pi * omega * g)
        / (2 * sp.pi * g) ** sp.Rational(2, 3)
    }
