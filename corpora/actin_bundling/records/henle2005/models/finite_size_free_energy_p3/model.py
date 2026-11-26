import math
from scipy.optimize import brentq, minimize_scalar
from scipy.special import k0e, k1e

M_TO_NM3 = 0.602214076


def _sigmoid(y):
    if y >= 0.0:
        return 1.0 / (1.0 + math.exp(-y))
    ey = math.exp(y)
    return ey / (1.0 + ey)


def _log1m_sigmoid(y):
    if y > 0.0:
        return -(y + math.log1p(math.exp(-y)))
    return -math.log1p(math.exp(y))


def _minimize_printed_free_energy_at(Ra, a_0, a, w, lam0, l_B, n_1_M, n_q_M, q, gamma, lam_star_frac, tol):
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

    def eval_at_logit(y):
        u = _sigmoid(y)
        lq = Lambda * u
        lam = lam0 - q * lq

        # Solve Eq. (7) for x = 2*l_B*lam_tot*Kratio given fixed lam
        def res_x(x):
            lam_st = alpha_st * (n1 * math.exp(x) - nm1 * math.exp(-x) + q * nq * math.exp(q * x))
            return 2.0 * l_B * (lam * M - lam_st) * Kratio - x

        x_lo, x_hi = -5.0, 20.0
        while res_x(x_hi) > 0.0:
            x_hi *= 2.0
        while res_x(x_lo) < 0.0:
            x_lo *= 2.0
        x = float(brentq(res_x, x_lo, x_hi, xtol=tol * 0.1, rtol=tol * 0.1))

        ls1 = n1 * alpha_st * math.exp(x)
        lsm1 = nm1 * alpha_st * math.exp(-x)
        lsq = nq * alpha_st * math.exp(q * x)
        lam_st = ls1 - lsm1 + q * lsq
        lam_tot = lam * M - lam_st

        log_1mu = _log1m_sigmoid(y)
        one_minus_u = math.exp(log_1mu)
        log_u = -math.log1p(math.exp(-y)) if y >= 0.0 else y - math.log1p(math.exp(y))
        log_lq_over_nqab = math.log(Lambda / (nq * alpha_b)) + log_u

        # Printed modified free energy per rod per unit length (Eqs. 2, 5 + p. 3 steric term)
        # Note exact algebraic cancellation of -lambda_q + lambda_q inside the bundle entropy sum
        F1_over_L = (
            (n1 + nm1 + nq) * alpha_st
            + ls1 * (x - 1.0)
            + lsm1 * (-x - 1.0)
            + lsq * (q * x - 1.0)
            + M * (lq * log_lq_over_nqab + (Lambda * one_minus_u) * log_1mu)
        )
        F2 = (
            3.6 * gamma / Ra
            + l_B * lam * lam * M / 4.0
            + l_B * lam_tot * lam_tot * Kratio / M
        )
        return F1_over_L / M + F2, lam_tot, lam

    # Directly minimize the printed modified free energy over logit y = ln(u/(1-u))
    # Coarse scan to bracket the minimum tightly before Brent bounded refinement
    y_est = max(0.0, 0.5 * l_B * q * (lam0 * lam_star_frac) * M)
    y_max = max(40.0, y_est + 60.0)
    grid = [-25.0 + (y_max + 25.0) * k / 40.0 for k in range(41)]
    vals = [eval_at_logit(yy)[0] for yy in grid]
    k_best = min(range(len(vals)), key=lambda k: vals[k])
    y_lo = grid[max(0, k_best - 1)]
    y_hi = grid[min(len(grid) - 1, k_best + 1)]
    res = minimize_scalar(lambda yy: eval_at_logit(yy)[0], bounds=(y_lo, y_hi), method="bounded", options={"xatol": tol})
    f_val, lam_tot, lam_rod = eval_at_logit(float(res.x))
    return f_val, lam_tot, lam_rod


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
    f_R, lam_tot, lam_rod = _minimize_printed_free_energy_at(inputs["R_over_a"], *args)
    f_ref, _, _ = _minimize_printed_free_energy_at(inputs["R_ref_over_a"], *args)
    return {
        "dF_over_NL": f_R - f_ref,
        "lambda_tot": lam_tot,
        "lambda_rod": lam_rod,
    }
