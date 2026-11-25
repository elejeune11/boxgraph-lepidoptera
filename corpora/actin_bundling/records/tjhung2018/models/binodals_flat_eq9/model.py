import numpy as np
from scipy.optimize import fsolve


def run(inputs, settings):
    A = float(inputs["A"])
    K = float(inputs["K"])
    zeta = float(inputs["zeta"])
    lam = float(inputs["lam"])
    n_quad = int(settings.get("n_quad", 64))
    xtol = float(settings.get("xtol", 1e-12))

    alpha = (zeta - 2.0 * lam) / K
    nodes, weights = np.polynomial.legendre.leggauss(n_quad)

    # Integrate dpsi/dphi = exp(alpha * phi) and dg/dphi = f'(phi) * exp(alpha * phi) on [0, phi]
    def psi_g(phi):
        half = 0.5 * phi
        xs = half * (nodes + 1.0)
        ex = np.exp(alpha * xs)
        psi_val = float(half * np.dot(weights, ex))
        g_val = float(half * np.dot(weights, A * (xs**3 - xs) * ex))
        return psi_val, g_val

    def res(x):
        p1, p2 = float(x[0]), float(x[1])
        mu1 = A * (p1**3 - p1)
        mu2 = A * (p2**3 - p2)
        psi1, g1 = psi_g(p1)
        psi2, g2 = psi_g(p2)
        return [mu1 - mu2, (psi1 * mu1 - g1) - (psi2 * mu2 - g2)]

    p1, p2 = fsolve(res, [-1.0, 1.0], xtol=xtol)
    p1, p2 = float(p1), float(p2)
    if p1 > p2:
        p1, p2 = p2, p1
    mu_bar = A * (p1**3 - p1)
    psi1, g1 = psi_g(p1)
    psi2, g2 = psi_g(p2)
    return {
        "phi_1": p1,
        "phi_2": p2,
        "mu_bar": mu_bar,
        "Delta_phi": p2 - p1,
        "Delta_psi": psi2 - psi1,
        "Delta_g": g2 - g1,
    }
