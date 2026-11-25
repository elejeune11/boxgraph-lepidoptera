import math
from scipy.optimize import fsolve


def _psi_g_printed(phi, A, K, zeta, lam):
    d = zeta - 2.0 * lam
    if abs(d) < 1e-12:
        if abs(zeta) < 1e-12:
            return phi, A * (-0.5 * phi**2 + 0.25 * phi**4)
        raise ValueError("Printed Eq. (B2) diverges at zeta - 2*lambda = 0 for zeta != 0")
    e = math.exp((d / K) * phi)
    psi = (K / d) * (e - 1.0)
    bracket = (
        -(d**3) * phi * (phi**2 - 1.0)
        + 6.0 * K**3
        - 6.0 * zeta * K**2 * phi
        + (d**2) * K * (3.0 * phi**2 - 1.0)
    )
    g = 6.0 * A * K**4 / d**4 - A * K**2 / d**2 - (A * K * e / d**4) * bracket
    return psi, g


def run(inputs, settings):
    A = float(inputs["A"])
    K = float(inputs["K"])
    combo = float(inputs["combo"])
    lam = float(inputs["lam"])
    zeta = combo + 2.0 * lam
    xtol = float(settings.get("xtol", 1e-12))

    def res(x):
        p1, p2 = float(x[0]), float(x[1])
        mu1 = A * (p1**3 - p1)
        mu2 = A * (p2**3 - p2)
        psi1, g1 = _psi_g_printed(p1, A, K, zeta, lam)
        psi2, g2 = _psi_g_printed(p2, A, K, zeta, lam)
        return [mu1 - mu2, (psi1 * mu1 - g1) - (psi2 * mu2 - g2)]

    p1, p2 = fsolve(res, [-1.0, 1.0], xtol=xtol)
    p1, p2 = float(p1), float(p2)
    if p1 > p2:
        p1, p2 = p2, p1
    return {"phi_1": p1, "phi_2": p2}
