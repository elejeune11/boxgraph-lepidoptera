import math
from scipy.optimize import fsolve


def _psi_g_corrected(phi, A, K, d):
    if abs(d) < 1e-6:
        # Smooth Taylor expansion around d = zeta - 2*lambda = 0
        a = d / K
        psi = phi + 0.5 * a * phi**2 + (a**2 / 6.0) * phi**3 + (a**3 / 24.0) * phi**4
        g = (
            A * (-0.5 * phi**2 + 0.25 * phi**4)
            + a * A * (-(1.0 / 3.0) * phi**3 + 0.2 * phi**5)
            + 0.5 * a**2 * A * (-0.25 * phi**4 + (1.0 / 6.0) * phi**6)
        )
        return psi, g
    e = math.exp((d / K) * phi)
    psi = (K / d) * (e - 1.0)
    bracket = (
        -(d**3) * phi * (phi**2 - 1.0)
        + 6.0 * K**3
        - 6.0 * d * K**2 * phi
        + (d**2) * K * (3.0 * phi**2 - 1.0)
    )
    g = 6.0 * A * K**4 / d**4 - A * K**2 / d**2 - (A * K * e / d**4) * bracket
    return psi, g


def run(inputs, settings):
    A = float(inputs["A"])
    K = float(inputs["K"])
    combo = float(inputs["combo"])
    xtol = float(settings.get("xtol", 1e-12))

    def res(x):
        p1, p2 = float(x[0]), float(x[1])
        mu1 = A * (p1**3 - p1)
        mu2 = A * (p2**3 - p2)
        psi1, g1 = _psi_g_corrected(p1, A, K, combo)
        psi2, g2 = _psi_g_corrected(p2, A, K, combo)
        return [mu1 - mu2, (psi1 * mu1 - g1) - (psi2 * mu2 - g2)]

    p1, p2 = fsolve(res, [-1.0, 1.0], xtol=xtol)
    p1, p2 = float(p1), float(p2)
    if p1 > p2:
        p1, p2 = p2, p1
    return {"phi_1": p1, "phi_2": p2}
