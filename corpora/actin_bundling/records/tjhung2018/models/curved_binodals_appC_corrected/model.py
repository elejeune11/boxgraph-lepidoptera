import math
import numpy as np
from scipy.optimize import fsolve


def _solve_appC_corrected(A, K, zeta, lam, d, R, n_grid, xtol):
    alpha = (zeta - 2.0 * lam) / K
    nodes, weights = np.polynomial.legendre.leggauss(64)

    def psi_g(phi):
        half = 0.5 * phi
        xs_q = half * (nodes + 1.0)
        ex = np.exp(alpha * xs_q)
        return float(half * np.dot(weights, ex)), float(half * np.dot(weights, A * (xs_q**3 - xs_q) * ex))

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

    xs = np.linspace(p1, p2, n_grid)
    w0 = np.empty_like(xs)
    for i, x in enumerate(xs):
        psix, gx = psi_g(float(x))
        w0[i] = max(0.0, -(2.0 / K) * ((psix * mu_bar - gx) - P1) * math.exp(-alpha * float(x)))
    sw0 = np.sqrt(w0)
    I0 = float(np.trapezoid(sw0, xs))
    S0 = math.exp(alpha * p2) * I0
    S1 = float(np.trapezoid(sw0 * np.exp(alpha * xs), xs))
    Delta1 = (d - 1.0) * zeta * I0 / R
    if abs(zeta - 2.0 * lam) < 1e-8:
        I_phi = float(np.trapezoid(xs * sw0, xs))
        sigma1 = K * I0 + zeta * (p2 * I0 - I_phi)
    else:
        sigma1 = (K / (zeta - 2.0 * lam)) * (zeta * S0 - 2.0 * lam * S1)

    def solve_pair(Delta, sigma, init):
        def res(x):
            pm, pp = float(x[0]), float(x[1])
            mu_m = A * (pm**3 - pm)
            mu_p = A * (pp**3 - pp) + Delta
            psim, gm = psi_g(pm)
            psip, gp = psi_g(pp)
            return [mu_m - mu_p, (psip * mu_m - gp) - (psim * mu_m - gm) - (d - 1.0) * sigma / R]
        pm, pp = fsolve(res, init, xtol=xtol)
        return float(pm), float(pp)

    pm1, pp1 = solve_pair(Delta1, sigma1, [p1, p2])
    muI1 = A * (pm1**3 - pm1)

    # Corrected O(1/R^2) step with phi' = -sqrt(w)
    xs2 = np.linspace(pm1, pp1, n_grid)
    sw_old = np.interp(xs2, xs, sw0, left=0.0, right=0.0)
    dx = xs2[1] - xs2[0]
    cum_sw = np.zeros_like(xs2)
    cum_sw[1:] = np.cumsum(0.5 * (sw_old[:-1] + sw_old[1:]) * dx)
    s_vals = (
        (2.0 / K) * (A * (xs2**3 - xs2) - muI1)
        + (2.0 * (d - 1.0) / R) * sw_old
        + (2.0 * (d - 1.0) * zeta / (K * R)) * cum_sw
    )
    integrand = np.exp(alpha * xs2) * s_vals
    cum_int = np.zeros_like(xs2)
    cum_int[1:] = np.cumsum(0.5 * (integrand[:-1] + integrand[1:]) * dx)
    w2 = np.maximum(0.0, np.exp(-alpha * xs2) * cum_int)
    sw2 = np.sqrt(w2)
    I0_2 = float(np.trapezoid(sw2, xs2))
    S0_2 = math.exp(alpha * pp1) * I0_2
    S1_2 = float(np.trapezoid(sw2 * np.exp(alpha * xs2), xs2))
    Delta2 = (d - 1.0) * zeta * I0_2 / R
    if abs(zeta - 2.0 * lam) < 1e-8:
        I_phi2 = float(np.trapezoid(xs2 * sw2, xs2))
        sigma2 = K * I0_2 + zeta * (pp1 * I0_2 - I_phi2)
    else:
        sigma2 = (K / (zeta - 2.0 * lam)) * (zeta * S0_2 - 2.0 * lam * S1_2)

    pm2, pp2 = solve_pair(Delta2, sigma2, [pm1, pp1])
    muI2 = A * (pm2**3 - pm2)
    return pm2, pp2, muI2, Delta2, sigma2


def run(inputs, settings):
    A = float(inputs["A"])
    K = float(inputs["K"])
    zeta = float(inputs["zeta"])
    lam = float(inputs["lam"])
    d = float(inputs["d"])
    R = float(inputs["R"])
    n_grid = int(settings.get("n_grid", 400))
    xtol = float(settings.get("xtol", 1e-12))

    pm, pp, muI, Delta, sigma = _solve_appC_corrected(A, K, zeta, lam, d, R, n_grid, xtol)
    dR = max(1e-3, 1e-3 * R)
    pm_hi, pp_hi, _, _, _ = _solve_appC_corrected(A, K, zeta, lam, d, R + dR, n_grid, xtol)
    pm_lo, pp_lo, _, _, _ = _solve_appC_corrected(A, K, zeta, lam, d, R - dR, n_grid, xtol)
    dpm_dR = (pm_hi - pm_lo) / (2.0 * dR)
    dpp_dR = (pp_hi - pp_lo) / (2.0 * dR)
    # Forward Ostwald: d(phi_-)/dR < 0 (+1); Reverse Ostwald: d(phi_-)/dR > 0 (-1)
    ostwald_sign = 1.0 if dpm_dR < 0.0 else -1.0
    return {
        "phi_minus": pm,
        "phi_plus": pp,
        "mu_I": muI,
        "Delta": Delta,
        "sigma": sigma,
        "dphi_minus_dR": dpm_dR,
        "dphi_plus_dR": dpp_dR,
        "ostwald_sign_from_slope": ostwald_sign,
    }
