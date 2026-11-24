import math
from scipy.optimize import brentq


def _mu(x):
    return x * (1.0 - x) * (1.0 - 2.0 * x)


def _f(x):
    return 0.5 * (x * (1.0 - x)) ** 2


def _solve_x_minus(mu_val, tol):
    if mu_val <= 0.0:
        return 1.0
    hi = 1.0 + max(1.0, (mu_val / 2.0) ** (1.0 / 3.0) + 1.0)
    while _mu(hi) < mu_val:
        hi *= 2.0
    xtol = max(tol * min(1.0, mu_val), 1e-15)
    return brentq(lambda xm: _mu(xm) - mu_val, 1.0, hi, xtol=xtol, rtol=tol)


def run(inputs, settings):
    b = inputs["b"]
    psi = inputs["psi"]
    gamma = inputs["gamma"]
    R = inputs["R"]
    tol = float(settings.get("tol", 1.0e-8))

    alpha = 2.0 * gamma / (b * psi**2 * R)
    if alpha <= 0.0:
        return {"phi_plus_B": 0.0, "phi_minus_B": psi}

    x_sp = (3.0 - math.sqrt(3.0)) / 6.0
    sub_tol = max(tol * 0.01, 1e-14)

    def residual(xp):
        mu_v = _mu(xp)
        xm = _solve_x_minus(mu_v, sub_tol)
        return (xm - xp) * mu_v + _f(xp) - _f(xm) - alpha

    if residual(x_sp) < 0.0:
        return {"phi_plus_B": None, "phi_minus_B": None}

    xtol = max(tol * min(1.0, alpha) * 0.01, 1e-15)
    xp = brentq(residual, 0.0, x_sp, xtol=xtol, rtol=tol)
    xm = _solve_x_minus(_mu(xp), sub_tol)
    return {
        "phi_plus_B": psi * xp,
        "phi_minus_B": psi * xm,
    }
