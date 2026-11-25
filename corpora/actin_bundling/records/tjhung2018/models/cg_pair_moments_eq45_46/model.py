import math
import numpy as np


def run(inputs, settings):
    d = float(inputs["d"])
    U_0 = float(inputs["U_0"])
    ell = float(inputs["ell"])
    n_quad = int(settings.get("n_quad", 64))

    # Surface area of unit sphere in R^d: Omega = 2 * pi^(d/2) / Gamma(d/2)
    Omega = 2.0 * (math.pi ** (0.5 * d)) / math.gamma(0.5 * d)
    # Direct Gauss-Legendre integration of Eqs. (45)-(46) for U(y) = U_0 exp(-y^2 / (2 ell^2))
    # over y in [0, 12*ell] (tail exp(-72) ~ 5e-32 < machine epsilon; integrand is entire).
    y_max = 12.0 * ell
    t_nodes, t_weights = np.polynomial.legendre.leggauss(n_quad)
    y = 0.5 * y_max * (t_nodes + 1.0)
    w = 0.5 * y_max * t_weights
    U_prime = -(U_0 / (ell * ell)) * y * np.exp(-0.5 * (y / ell) ** 2)
    int_B = float(np.dot(w, (y ** (d + 2.0)) * U_prime))
    int_C = float(np.dot(w, (y ** d) * U_prime))
    B = -(Omega / (2.0 * d * (d + 2.0))) * int_B
    C = -(Omega / d) * int_C
    return {"B": B, "C": C}


def symbolic():
    import sympy as sp

    d, U_0, ell = sp.symbols("d U_0 ell", positive=True)
    # Exact analytical evaluation of Eqs. (45)-(46) for U(y) = U_0 * exp(-y^2 / (2 ell^2)):
    # C = int U(r) d^d r = U_0 * (2 pi ell^2)^(d/2)
    # B = (1 / (2 d)) int r^2 U(r) d^d r = (1/2) * ell^2 * C
    C = U_0 * (2 * sp.pi * ell**2) ** (d / 2)
    B = sp.Rational(1, 2) * ell**2 * C
    return {"B": B, "C": C}
