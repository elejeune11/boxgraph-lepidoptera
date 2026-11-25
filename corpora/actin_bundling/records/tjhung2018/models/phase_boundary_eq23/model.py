import math
import numpy as np
from scipy.optimize import brentq, fsolve


def _sigma_droplet(A, K, zeta, lam, n_quad):
    alpha = (zeta - 2.0 * lam) / K
    nodes, weights = np.polynomial.legendre.leggauss(n_quad)

    def psi_g(phi):
        half = 0.5 * phi
        xs = half * (nodes + 1.0)
        ex = np.exp(alpha * xs)
        return float(half * np.dot(weights, ex)), float(half * np.dot(weights, A * (xs**3 - xs) * ex))

    def res(x):
        p1, p2 = float(x[0]), float(x[1])
        mu1 = A * (p1**3 - p1)
        mu2 = A * (p2**3 - p2)
        psi1, g1 = psi_g(p1)
        psi2, g2 = psi_g(p2)
        return [mu1 - mu2, (psi1 * mu1 - g1) - (psi2 * mu2 - g2)]

    p1, p2 = fsolve(res, [-1.0, 1.0], xtol=1e-12)
    p1, p2 = float(p1), float(p2)
    if p1 > p2:
        p1, p2 = p2, p1
    mu_bar = A * (p1**3 - p1)
    psi1, g1 = psi_g(p1)
    P1 = psi1 * mu_bar - g1

    mid = 0.5 * (p1 + p2)
    half = 0.5 * (p2 - p1)
    xs = mid + half * nodes
    w0 = np.empty_like(xs)
    for i, x in enumerate(xs):
        psix, gx = psi_g(float(x))
        val = -(2.0 / K) * ((psix * mu_bar - gx) - P1) * math.exp(-alpha * float(x))
        w0[i] = max(0.0, val)
    sw0 = np.sqrt(w0)
    I0 = float(half * np.dot(weights, sw0))
    S0 = math.exp(alpha * p2) * I0
    S1 = float(half * np.dot(weights, sw0 * np.exp(alpha * xs)))
    if abs(zeta - 2.0 * lam) < 1e-8:
        I_phi = float(half * np.dot(weights, xs * sw0))
        return K * I0 + zeta * (p2 * I0 - I_phi)
    return (K / (zeta - 2.0 * lam)) * (zeta * S0 - 2.0 * lam * S1)


def run(inputs, settings):
    A = float(inputs["A"])
    K = float(inputs["K"])
    zeta = float(inputs["zeta"])
    n_quad = int(settings.get("n_quad", 128))
    xtol = float(settings.get("xtol", 1e-10))

    if zeta == 0.0:
        raise ValueError("No reverse-Ostwald boundary exists at zeta = 0 (AMB lies strictly in Region A)")
    sign = 1.0 if zeta < 0.0 else -1.0
    z_eff = -abs(zeta)
    # For droplets (z_eff < 0): sigma(lam) > 0 for large positive lam, < 0 for large negative lam
    lo, hi = -8.0 * K, 4.0 * K
    while _sigma_droplet(A, K, z_eff, lo, n_quad) > 0.0:
        lo *= 2.0
    lam_drop = brentq(lambda l: _sigma_droplet(A, K, z_eff, l, n_quad), lo, hi, xtol=xtol)
    # If input zeta > 0, return Region A/C bubble boundary lambda*_bubble = -lam_drop
    return {"lam_star": sign * float(lam_drop)}
