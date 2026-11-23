import math
from scipy.optimize import brentq


def _chi_c2(phi):
    s = math.sqrt(phi)
    return (2.0 * s**3 + 3.0 * phi - 1.0) / (9.0 * (s + 2.0))


def run(inputs, settings):
    phi_c = brentq(_chi_c2, 1e-9, 10.0)
    return {"chi_c2": _chi_c2(inputs["phi"]), "phi_tricritical": phi_c}


def symbolic():
    import sympy as sp

    phi = sp.Symbol("phi", positive=True)
    return {
        "chi_c2": (2 * phi ** sp.Rational(3, 2) + 3 * phi - 1) / (9 * (sp.sqrt(phi) + 2)),
        "phi_tricritical": sp.Rational(1, 4),
    }
