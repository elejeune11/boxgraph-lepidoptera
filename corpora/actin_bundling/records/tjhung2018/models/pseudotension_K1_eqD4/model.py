import math
import numpy as np
from scipy.optimize import fsolve


def run(inputs, settings):
    A = float(inputs["A"])
    K = float(inputs["K"])
    K_1 = float(inputs["K_1"])
    zeta = float(inputs["zeta"])
    lam = float(inputs["lam"])
    n_quad = int(settings.get("n_quad", 128))
    xtol = float(settings.get("xtol", 1e-12))

    nu = lam - K_1 - 0.5 * zeta
    nodes, weights = np.polynomial.legendre.leggauss(n_quad)

    def psi_bar_arr(x):
        if abs(K_1) < 1e-10:
            return np.exp(-2.0 * nu * x / K)
        return np.power(1.0 + 2.0 * K_1 * x / K, -nu / K_1)

    def psi_prime_arr(x):
        # d(psi_{K_1})/dphi = (K / K(phi)) * psi_bar(phi)
        return (K / (K + 2.0 * K_1 * x)) * psi_bar_arr(x)

    def psi_g(phi):
        half = 0.5 * phi
        xs = half * (nodes + 1.0)
        pp = psi_prime_arr(xs)
        return float(half * np.dot(weights, pp)), float(half * np.dot(weights, A * (xs**3 - xs) * pp))

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
    P1 = psi1 * mu_bar - g1

    # From 1D stationary Eq. (D1) multiplied by d_x psi_{K_1}:
    # d_x [-(1/2) K(phi) psi'_{K_1}(phi) phi'^2] = d_x [P(phi)] -> (1/2) K * psi_bar(phi) * w0(phi) = -(P(phi) - P(phi_1))
    mid = 0.5 * (p1 + p2)
    half = 0.5 * (p2 - p1)
    xs = mid + half * nodes
    pb = psi_bar_arr(xs)
    w0 = np.empty_like(xs)
    for i, x in enumerate(xs):
        psix, gx = psi_g(float(x))
        Px = psix * mu_bar - gx
        w0[i] = max(0.0, -(2.0 / (K * float(pb[i]))) * (Px - P1))
    sw0 = np.sqrt(w0)
    I0 = float(half * np.dot(weights, sw0))
    pb2 = float(psi_bar_arr(np.array([p2]))[0])
    S2 = pb2 * I0
    S3 = float(half * np.dot(weights, pb * sw0))

    if abs(nu) < 1e-8:
        # Smooth nu -> 0 limit (zeta = 2*lam - 2*K_1)
        if abs(K_1) < 1e-10:
            I_log = float(half * np.dot(weights, xs * sw0))
            log2 = p2
        else:
            I_log = float(half * np.dot(weights, (K / (2.0 * K_1)) * np.log(1.0 + 2.0 * K_1 * xs / K) * sw0))
            log2 = (K / (2.0 * K_1)) * math.log(1.0 + 2.0 * K_1 * p2 / K)
        sigma_K1 = K * I0 + zeta * (log2 * I0 - I_log)
    else:
        sigma_K1 = -(K / (2.0 * nu)) * (zeta * S2 - (2.0 * lam - 2.0 * K_1) * S3)

    return {
        "phi_1": p1,
        "phi_2": p2,
        "S_2": S2,
        "S_3": S3,
        "sigma_K1": sigma_K1,
    }
