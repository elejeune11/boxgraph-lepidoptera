import math
from scipy.optimize import brentq
from scipy.special import k0e, k1e

M_TO_NM3 = 0.602214076


def _sigmoid(y):
    if y >= 0.0:
        return 1.0 / (1.0 + math.exp(-y))
    ey = math.exp(y)
    return ey / (1.0 + ey)


def _log1m_sigmoid(y):
    # ln(1 - sigmoid(y)) = -ln(1 + exp(y))
    if y > 0.0:
        return -(y + math.log1p(math.exp(-y)))
    return -math.log1p(math.exp(y))


def _solve_finite_at(Ra, a_0, a, w, lam0, l_B, n_1_M, n_q_M, q, gamma, lam_star_frac, tol):
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

    Lambda = lam0 * (1.0 - lam_star_frac) / q
    log_L_over_nqab = math.log(Lambda / (nq * alpha_b))

    def lam_q_of_x(x):
        c = l_B * q * M / 2.0
        # Corrected LHS Eq. (6) in logit variable y = ln(u/(1-u)), u = lam_q/Lambda:
        # y + c*q*Lambda*sigmoid(y) = c*lam0 + q*x - ln(Lambda/(n_q*alpha_b))
        rhs = c * lam0 + q * x - log_L_over_nqab
        y_lo = rhs - c * q * Lambda - 5.0
        y_hi = rhs + 5.0
        y = float(brentq(lambda yy: yy + c * q * Lambda * _sigmoid(yy) - rhs, y_lo, y_hi, xtol=tol, rtol=tol))
        return Lambda * _sigmoid(y), y

    def residual(x):
        lq, _ = lam_q_of_x(x)
        lam = lam0 - q * lq
        lam_st = alpha_st * (n1 * math.exp(x) - nm1 * math.exp(-x) + q * nq * math.exp(q * x))
        lam_tot = lam * M - lam_st
        return 2.0 * l_B * lam_tot * Kratio - x

    x_lo, x_hi = -5.0, 20.0
    while residual(x_hi) > 0.0:
        x_hi *= 2.0
    while residual(x_lo) < 0.0:
        x_lo *= 2.0
    x = float(brentq(residual, x_lo, x_hi, xtol=tol, rtol=tol))

    lq, y = lam_q_of_x(x)
    log_1mu = _log1m_sigmoid(y)
    one_minus_u = math.exp(log_1mu)
    lam = lam0 - q * lq

    ls1 = n1 * alpha_st * math.exp(x)
    lsm1 = nm1 * alpha_st * math.exp(-x)
    lsq = nq * alpha_st * math.exp(q * x)
    lam_st = ls1 - lsm1 + q * lsq
    lam_tot = lam * M - lam_st

    F1_over_L = (
        (n1 + nm1 + nq) * alpha_st
        + ls1 * (x - 1.0)
        + lsm1 * (-x - 1.0)
        + lsq * (q * x - 1.0)
        + lq * M * (math.log(lq / (nq * alpha_b)) - 1.0)
        + M * ((Lambda * one_minus_u) * log_1mu + lq)
    )
    F2 = (
        3.6 * gamma / Ra
        + l_B * lam * lam * M / 4.0
        + l_B * lam_tot * lam_tot * Kratio / M
    )
    f_over_NL = F1_over_L / M + F2
    return f_over_NL, lam_tot, lam


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
        inputs["lambda_star_frac"],
        tol,
    )
    f_R, lam_tot, lam_rod = _solve_finite_at(inputs["R_over_a"], *args)
    f_ref, _, _ = _solve_finite_at(inputs["R_ref_over_a"], *args)
    return {
        "dF_over_NL": f_R - f_ref,
        "lambda_tot": lam_tot,
        "lambda_rod": lam_rod,
    }
