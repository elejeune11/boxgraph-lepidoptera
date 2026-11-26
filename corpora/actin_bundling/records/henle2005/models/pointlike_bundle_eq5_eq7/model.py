import math
from scipy.optimize import brentq
from scipy.special import k0e, k1e

M_TO_NM3 = 0.602214076


def _solve_pointlike_at(Ra, a_0, a, w, lam0, l_B, n_1_M, n_q_M, q, gamma, tol):
    R = Ra * a
    M = Ra * Ra
    alpha_b = math.pi * (a * a - a_0 * a_0)
    alpha_st = 2.0 * math.pi * R * w
    n1 = n_1_M * M_TO_NM3
    nq = n_q_M * M_TO_NM3
    nm1 = n1 + q * nq
    kappa = math.sqrt(4.0 * math.pi * l_B * (n1 + nm1 + q * q * nq))

    kR = kappa * R
    # Use exponentially scaled Bessel functions k0e/k1e = k0/k1 to avoid overflow/underflow at large kappa*R
    Kratio = float(k0e(kR) / (kR * k1e(kR)))

    def u_of_x(x):
        # Solve Eq. (6) directly for u = lam * M in (-inf, lam0 * M) (O(10) at all R, even R -> inf):
        # ln((lam0 - u/M)/(q * nq * alpha_b)) - (l_B * q / 2) * u - q * x = 0
        c_u = l_B * q / 2.0
        log_denom = math.log(q * nq * alpha_b)
        u_cap = lam0 * M

        def f(u):
            return math.log(lam0 - u / M) - log_denom - c_u * u - q * x

        # Asymptotic scale: u ~ (ln(lam0/(q*nq*alpha_b)) - q*x) / c_u
        u_hi = min(u_cap * (1.0 - 1e-15), max(10.0, (math.log(lam0) - log_denom - q * x) / c_u + 10.0))
        while f(u_hi) > 0.0:
            if u_hi < u_cap * 0.5:
                u_hi = min(u_cap * (1.0 - 1e-15), u_hi * 2.0 + 10.0)
            else:
                u_hi = u_cap - 0.5 * (u_cap - u_hi)
                if u_cap - u_hi < 1e-300 * u_cap:
                    break
        u_lo = 0.0
        if f(u_lo) < 0.0:
            u_lo = -10.0
            while f(u_lo) < 0.0:
                u_lo *= 2.0
        u = float(brentq(f, u_lo, u_hi, xtol=tol, rtol=tol))
        lam = u / M
        lq = (lam0 - lam) / q
        return u, lam, lq

    def lam_st_of_x(x):
        # Exact cancellation-free Stern charge using bulk neutrality n_{-1} = n_1 + q*n_q:
        # n1*(e^x - e^-x) + q*nq*(e^(q*x) - e^-x) = 2*n1*sinh(x) + q*nq*(expm1(q*x) - expm1(-x))
        return alpha_st * (2.0 * n1 * math.sinh(x) + q * nq * (math.expm1(q * x) - math.expm1(-x)))

    def residual_ltot(ltot):
        x = 2.0 * l_B * ltot * Kratio
        u, _, _ = u_of_x(x)
        return (u - lam_st_of_x(x)) - ltot

    lt_lo, lt_hi = -1.0, 20.0
    while residual_ltot(lt_hi) > 0.0:
        lt_hi *= 2.0
    while residual_ltot(lt_lo) < 0.0:
        lt_lo *= 2.0
    lam_tot = float(brentq(residual_ltot, lt_lo, lt_hi, xtol=tol, rtol=tol))

    x = 2.0 * l_B * lam_tot * Kratio
    u, lam, lq = u_of_x(x)
    ls1 = n1 * alpha_st * math.exp(x)
    lsm1 = nm1 * alpha_st * math.exp(-x)
    lsq = nq * alpha_st * math.exp(q * x)

    # Free energy per rod per unit length: (a^2/R^2)*(F_1/L) + F_2 (Eqs. 2 & 5)
    # Use cancellation-free Stern entropy sum_s [n_s*alpha_st*(1 + e^{s*x}*(s*x - 1))]
    def _stern_term(ns, sx):
        if abs(sx) < 1e-4:
            # 1 + e^y*(y - 1) = y^2/2 + y^3/3 + y^4/8 + ...
            return ns * alpha_st * (0.5 * sx * sx + (sx ** 3) / 3.0 + (sx ** 4) / 8.0)
        return ns * alpha_st * (1.0 + math.exp(sx) * (sx - 1.0))

    F1_over_L = (
        _stern_term(n1, x)
        + _stern_term(nm1, -x)
        + _stern_term(nq, q * x)
        + lq * M * (math.log(lq / (nq * alpha_b)) - 1.0)
    )
    F2 = (
        3.6 * gamma / Ra
        + l_B * u * u / (4.0 * M)
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
        tol,
    )
    f_R, lam_tot, lam_rod = _solve_pointlike_at(inputs["R_over_a"], *args)
    f_ref, _, _ = _solve_pointlike_at(inputs["R_ref_over_a"], *args)
    return {
        "dF_over_NL": f_R - f_ref,
        "lambda_tot": lam_tot,
        "lambda_rod": lam_rod,
    }
