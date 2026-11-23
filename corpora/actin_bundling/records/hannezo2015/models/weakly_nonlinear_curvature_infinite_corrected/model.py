import math
from scipy.optimize import brentq


def _chi_c2(phi):
    s = math.sqrt(phi)
    return (1.0 + s) * ((s + phi) ** 2) * (2.0 * s - 1.0) / (9.0 * (s ** 3))


def run(inputs, settings):
    phi_c = brentq(_chi_c2, 1e-9, 10.0)
    return {"chi_c2": _chi_c2(inputs["phi"]), "phi_tricritical": phi_c}


def symbolic():
    import sympy as sp

    phi = sp.Symbol("phi", positive=True)
    s = sp.sqrt(phi)
    return {
        "chi_c2": (1 + s) * (s + phi) ** 2 * (2 * s - 1) / (9 * phi ** sp.Rational(3, 2)),
        "phi_tricritical": sp.Rational(1, 4),
    }
