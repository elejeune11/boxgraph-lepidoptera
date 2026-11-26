import math


def run(inputs, settings):
    b = inputs["b"]
    phi_max = inputs["phi_max"]
    xi = inputs["xi"]
    # Algebraically identical to 0.5*(b + phi_max) - 0.5*|phi_max - b|*sqrt(1 + xi/(b - phi_max)^2),
    # written in conjugate form to avoid catastrophic cancellation at large |phi_max - b|:
    d = phi_max - b
    s = math.sqrt(d * d + xi)
    if d >= 0.0:
        phi = b - 0.5 * xi / (d + s)
    else:
        phi = phi_max - 0.5 * xi / (-d + s)
    return {"phi": phi}


def symbolic():
    import sympy as sp

    b, phi_max, xi = sp.symbols("b phi_max xi", real=True)
    return {"phi": sp.Rational(1, 2) * (b + phi_max) - sp.Rational(1, 2) * sp.sqrt((phi_max - b) ** 2 + xi)}
