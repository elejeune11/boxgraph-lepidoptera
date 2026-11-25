import math
import numpy as np
from scipy.optimize import fsolve


def run(inputs, settings):
    A = float(inputs["A"])
    K = float(inputs["K"])
    zeta = float(inputs["zeta"])
    lam = float(inputs["lam"])
    d = float(inputs["d"])
    R = float(inputs["R"])
    n_quad = int(settings.get("n_quad", 128))
    xtol = float(settings.get("xtol", 1e-12))

    alpha = (zeta - 2.0 * lam) / K
    nodes, weights = np.polynomial.legendre.leggauss(n_quad)

    def psi_g(phi):
        half = 0.5 * phi
        xs = half * (nodes + 1.0)
        ex = np.exp(alpha * xs)
        return float(half * np.dot(weights, ex)), float(half * np.dot(weights, A * (xs**3 - xs) * ex))

    def res_flat(x):
        p1, p2 = float(x[0]), float(x[1])
        mu1 = A * (p1**3 - p1)
        mu2 = A * (p2**3 - p2)
        psi1, g1 = psi_g(p1)
        psi2, g2 = psi_g(p2)
        return [mu1 - mu2, (psi1 * mu1 - g1) - (psi2 * mu2 - g2)]

    p1, p2 = fsolve(res_flat, [-1.0, 1.0], xtol=xtol)
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
    Delta = (d - 1.0) * zeta * I0 / R
    if abs(zeta - 2.0 * lam) < 1e-8:
        I_phi = float(half * np.dot(weights, xs * sw0))
        sigma = K * I0 + zeta * (p2 * I0 - I_phi)
    else:
        sigma = (K / (zeta - 2.0 * lam)) * (zeta * S0 - 2.0 * lam * S1)

    def res_curved(x):
        pm, pp = float(x[0]), float(x[1])
        mu_m = A * (pm**3 - pm)
        mu_p = A * (pp**3 - pp) + Delta
        psim, gm = psi_g(pm)
        psip, gp = psi_g(pp)
        Pp = psip * mu_m - gp
        Pm = psim * mu_m - gm
        return [mu_m - mu_p, Pp - Pm - (d - 1.0) * sigma / R]

    pm, pp = fsolve(res_curved, [p1, p2], xtol=xtol)
    pm, pp = float(pm), float(pp)
    mu_I = A * (pm**3 - pm)
    return {
        "phi_minus": pm,
        "phi_plus": pp,
        "mu_I": mu_I,
        "Delta": Delta,
        "sigma": sigma,
    }
