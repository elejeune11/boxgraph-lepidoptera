import math
from scipy.optimize import brentq
from scipy.special import k0e, k1e

M_TO_NM3 = 0.602214076


def _phi_frustrated(b, phi_max, xi):
    d = phi_max - b
    s = math.sqrt(d * d + xi)
    if d >= 0.0:
        return b - 0.5 * xi / (d + s)
    return phi_max - 0.5 * xi / (-d + s)


def _solve_frustrated_at(Ra, a_0, a, w, lam0, l_B, n_1_M, n_q_M, q, gamma, phi_max, xi, tol):
    R = Ra * a
    M = Ra * Ra
    alpha_b = math.pi * (a * a - a_0 * a_0)
    alpha_st = 2.0 * math.pi * R * w
    n1 = n_1_M * M_TO_NM3
    nq = n_q_M * M_TO_NM3
    nm1 = n1 + q * nq
    kappa = math.sqrt(4.0 * math.pi * l_B * (n1 + nm1 + q * q * nq))

    kR = kappa * R
    Kratio = float(k0e(kR) / (kR * k1e(kR)))

    def lam_of_x(x):
        c = l_B * q * M / 2.0
        log_denom = math.log(q * nq * alpha_b)

        def f(lam):
            return math.log(lam0 - lam) - log_denom - c * lam - q * x

        hi = lam0 * (1.0 - 1e-15)
        while f(hi) > 0.0:
            hi = lam0 - 0.5 * (lam0 - hi)
            if lam0 - hi < 1e-300:
                break
        lo = 0.0
        if f(lo) < 0.0:
            lo = -1.0
            while f(lo) < 0.0:
                lo *= 2.0
        lam = float(brentq(f, lo, hi, xtol=tol, rtol=tol))
        lq = (lam0 - lam) / q
        return lam, lq

    def residual(x):
        lam, _ = lam_of_x(x)
        lam_st = alpha_st * (n1 * math.exp(x) - nm1 * math.exp(-x) + q * nq * math.exp(q * x))
        return 2.0 * l_B * (lam * M - lam_st) * Kratio - x

    x_lo, x_hi = -5.0, 20.0
    while residual(x_hi) > 0.0:
        x_hi *= 2.0
    while residual(x_lo) < 0.0:
        x_lo *= 2.0
    x = float(brentq(residual, x_lo, x_hi, xtol=tol, rtol=tol))

    lam, lq = lam_of_x(x)
    ls1 = n1 * alpha_st * math.exp(x)
    lsm1 = nm1 * alpha_st * math.exp(-x)
    lsq = nq * alpha_st * math.exp(q * x)
    lam_tot = lam * M - (ls1 - lsm1 + q * lsq)

    F1_over_L = (
        (n1 + nm1 + nq) * alpha_st
        + ls1 * (x - 1.0)
        + lsm1 * (-x - 1.0)
        + lsq * (q * x - 1.0)
        + lq * M * (math.log(lq / (nq * alpha_b)) - 1.0)
    )
    b = 3.0 - 3.6 / Ra
    F_attr_over_ML = -gamma * _phi_frustrated(b, phi_max, xi)
    F_es_over_ML = l_B * lam * lam * M / 4.0 + l_B * lam_tot * lam_tot * Kratio / M
    return F1_over_L / M + F_attr_over_ML + F_es_over_ML


def run(inputs, settings):
    tol = float(settings.get("tol", 1.0e-10))
    lam0 = 1.0 / inputs["inv_lambda_0"]
    args = (
        inputs["a_0"],
        inputs["a"],
        inputs["w"],
        lam0,
        inputs["l_B"],
        inputs["n_1"],
        inputs["n_q"],
        inputs["q"],
        inputs["gamma"],
        inputs["phi_max"],
        inputs["xi"],
        tol,
    )
    f_R = _solve_frustrated_at(inputs["R_over_a"], *args)
    f_ref = _solve_frustrated_at(inputs["R_ref_over_a"], *args)
    return {"dF_over_NL": f_R - f_ref}
