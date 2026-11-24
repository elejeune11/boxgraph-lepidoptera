import math


def run(inputs, settings):
    R = inputs["R"]
    a = inputs["alpha"]
    s = math.sin(a)
    c = math.cos(a)
    f = s - a * c
    g = f - (1.0 / 3.0) * (s**3)
    L = 2.0 * R * s
    r = 2.0 * R * (1.0 - c)
    return {
        "L": L,
        "r": r,
        "aspect_ratio": s / (1.0 - c),
        "A": 4.0 * math.pi * (R**2) * f,
        "A_omega": 2.0 * math.pi * (R**2) * g,
        "V": 2.0 * math.pi * (R**3) * g,
    }


def symbolic():
    import sympy as sp

    R, alpha = sp.symbols("R alpha", positive=True)
    f = sp.sin(alpha) - alpha * sp.cos(alpha)
    g = f - sp.Rational(1, 3) * sp.sin(alpha) ** 3
    return {
        "L": 2 * R * sp.sin(alpha),
        "r": 2 * R * (1 - sp.cos(alpha)),
        "aspect_ratio": sp.sin(alpha) / (1 - sp.cos(alpha)),
        "A": 4 * sp.pi * R**2 * f,
        "A_omega": 2 * sp.pi * R**2 * g,
        "V": 2 * sp.pi * R**3 * g,
    }
