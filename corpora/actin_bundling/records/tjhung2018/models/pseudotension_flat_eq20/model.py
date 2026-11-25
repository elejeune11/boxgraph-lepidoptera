import math
import numpy as np
from scipy.optimize import fsolve


def _droplet_pseudotension(A, K, zeta, lam, n_quad, xtol):
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

    p1, p2 = fsolve(res, [-1.0, 1.0], xtol=xtol)
    p1, p2 = float(p1), float(p2)
    if p1 > p2:
        p1, p2 = p2, p1
    mu_bar = A * (p1**3 - p1)
    psi1, g1 = psi_g(p1)
    P1 = psi1 * mu_bar - g1

    # Flat-interface w_0(x) = (phi')^2 = -(2/K) exp(-alpha*x) [P(x) - P(phi_1)] on [phi_1, phi_2]
    # Gauss-Legendre quadrature on [p1, p2]
    mid = 0.5 * (p1 + p2)
    half = 0.5 * (p2 - p1)
    xs = mid + half * nodes
    w0 = np.empty_like(xs)
    for i, x in enumerate(xs):
        psix, gx = psi_g(float(x))
        Px = psix * mu_bar - gx
        val = -(2.0 / K) * (Px - P1) * math.exp(-alpha * float(x))
        w0[i] = max(0.0, val)
    sw0 = np.sqrt(w0)
    I0 = float(half * np.dot(weights, sw0))
    S0 = math.exp(alpha * p2) * I0
    S1 = float(half * np.dot(weights, sw0 * np.exp(alpha * xs)))

    if abs(zeta - 2.0 * lam) < 1e-8:
        # Smooth analytical limit of Eq. (20) as zeta - 2*lambda -> 0:
        # sigma = K * I0 + zeta * (phi_2 * I0 - int phi * sqrt(w0) dphi)
        I_phi = float(half * np.dot(weights, xs * sw0))
        sigma = K * I0 + zeta * (p2 * I0 - I_phi)
    else:
        sigma = (K / (zeta - 2.0 * lam)) * (zeta * S0 - 2.0 * lam * S1)
    margin = -zeta * abs(1.0 - S0 / S1) - abs(zeta - 2.0 * lam)
    return S0, S1, sigma, margin


def run(inputs, settings):
    A = float(inputs["A"])
    K = float(inputs["K"])
    zeta = float(inputs["zeta"])
    lam = float(inputs["lam"])
    n_quad = int(settings.get("n_quad", 128))
    xtol = float(settings.get("xtol", 1e-12))

    S0, S1, sigma, margin = _droplet_pseudotension(A, K, zeta, lam, n_quad, xtol)
    # By duality (phi, lambda, zeta) -> -(phi, lambda, zeta), bubble pseudotension equals
    # droplet pseudotension evaluated at (-zeta, -lambda).
    _, _, sigma_bubble, _ = _droplet_pseudotension(A, K, -zeta, -lam, n_quad, xtol)
    return {
        "S_0": S0,
        "S_1": S1,
        "sigma": sigma,
        "sigma_bubble": sigma_bubble,
        "reverse_margin": margin,
        "ostwald_sign_droplet": 1.0 if sigma > 0.0 else -1.0,
        "ostwald_sign_bubble": 1.0 if sigma_bubble > 0.0 else -1.0,
    }
