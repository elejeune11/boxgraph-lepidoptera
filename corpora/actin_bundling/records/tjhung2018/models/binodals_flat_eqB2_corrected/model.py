import math
from scipy.optimize import fsolve


def _psi_g_exact(phi, A, alpha):
    if abs(alpha) < 1e-3:
        psi_val = 0.0
        g_val = 0.0
        fact = 1.0
        apow = 1.0
        for n in range(14):
            if n > 0:
                fact *= n
                apow *= alpha
            psi_val += apow * (phi ** (n + 1)) / (fact * (n + 1))
            g_val += A * (apow / fact) * ((phi ** (n + 4)) / (n + 4) - (phi ** (n + 2)) / (n + 2))
        return psi_val, g_val
    e = math.exp(alpha * phi)
    psi_val = (e - 1.0) / alpha
    bracket = (
        -(alpha**3) * phi * (phi**2 - 1.0)
        + 6.0
        - 6.0 * alpha * phi
        + (alpha**2) * (3.0 * phi**2 - 1.0)
    )
    g_val = A * (6.0 / alpha**4 - 1.0 / alpha**2 - (e / alpha**4) * bracket)
    return psi_val, g_val


def run(inputs, settings):
    A = float(inputs["A"])
    K = float(inputs["K"])
    zeta = float(inputs["zeta"])
    lam = float(inputs["lam"])
    xtol = float(settings.get("xtol", 1e-12))
    alpha = (zeta - 2.0 * lam) / K

    def res(x):
        p1, p2 = float(x[0]), float(x[1])
        mu1 = A * (p1**3 - p1)
        mu2 = A * (p2**3 - p2)
        psi1, g1 = _psi_g_exact(p1, A, alpha)
        psi2, g2 = _psi_g_exact(p2, A, alpha)
        return [mu1 - mu2, (psi1 * mu1 - g1) - (psi2 * mu2 - g2)]

    p1, p2 = fsolve(res, [-1.0, 1.0], xtol=xtol)
    p1, p2 = float(p1), float(p2)
    if p1 > p2:
        p1, p2 = p2, p1
    mu_bar = A * (p1**3 - p1)
    return {"phi_1": p1, "phi_2": p2, "mu_bar": mu_bar}
