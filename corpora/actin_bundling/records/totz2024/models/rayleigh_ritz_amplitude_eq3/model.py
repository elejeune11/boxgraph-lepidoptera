import math
from scipy.optimize import brentq


def run(inputs, settings):
    g = inputs["g"]
    p = inputs["p"]
    Y = inputs["Y"]
    nu = inputs["nu"]
    h = inputs["h"]
    L = inputs["L"]
    tol = float(settings.get("tol", 1.0e-6))

    g_c = (math.pi ** 2) / (3.0 * (1.0 - nu ** 2)) * ((h / L) ** 2)
    if p == 0.0:
        if g <= g_c:
            return {"a": 0.0}
        return {"a": (2.0 / (math.sqrt(3.0) * h)) * math.sqrt(g / g_c - 1.0)}

    # Printed Eq. 3 / E49 in SI base units (m, Pa):
    # f(a) = (3*h^2/4)*a^3 + (1 - g/g_c)*a - p*L^2/(pi^2*Y*g_c*h^3) = 0
    c3 = 0.75 * (h ** 2)
    c1 = 1.0 - g / g_c
    c0 = -p * (L ** 2) / ((math.pi ** 2) * Y * g_c * (h ** 3))

    def f(a):
        return c3 * (a ** 3) + c1 * a + c0

    lo = 1.0e-18
    while f(lo) >= 0.0:
        lo *= 0.1
    hi = 1.0e-9
    while f(hi) <= 0.0:
        hi *= 10.0

    a_sol = brentq(f, lo, hi, rtol=tol, xtol=1.0e-30, maxiter=500)
    return {"a": float(a_sol)}
