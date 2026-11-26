import math
from scipy.optimize import brentq

_EULER_GAMMA = 0.5772156649015329


def _log1pexp(x):
    if x > 50.0:
        return x
    if x < -50.0:
        return math.exp(x)
    return math.log1p(math.exp(x))


def run(inputs, settings):
    D = inputs["D"]
    L = inputs["L"]
    delta = inputs["delta"]
    kappa_inv = inputs["kappa_inv"]
    n_perp = inputs["n_perp"]
    Gamma_par = inputs["Gamma_par"]
    Gamma_perp = inputs["Gamma_perp"]
    eps_1 = inputs["eps_1"]
    eps_par = inputs["eps_par"]

    D_eff = D + kappa_inv * (_EULER_GAMMA + math.log(Gamma_perp))
    n_par = L / delta
    t1 = 8.0 * (D ** 2 + (kappa_inv ** 2) * (math.log(Gamma_par * L / D) ** 2)) / (D_eff * L)
    c2 = 16.0 * (D * delta / (D_eff * L))

    def B0_corr(mu_l):
        u_par_tot = -n_par * _log1pexp(mu_l - 2.0 * eps_par) + 2.0 * n_par * _log1pexp(mu_l - eps_1)
        Y = Gamma_par * (L / D) - abs(u_par_tot)
        if abs(Y) < 1e-8:
            br = 0.5 * Y
        elif Y < -500.0:
            return -1e200
        else:
            br = 1.0 + (math.exp(-Y) - 1.0) / Y
        return t1 + c2 * br

    mu_lo = 2.0 * eps_par - 20.0
    mu_hi = 2.0 * eps_par + 5.0
    grid = [mu_lo + (mu_hi - mu_lo) * i / 400.0 for i in range(401)]
    vals = [B0_corr(m) for m in grid]
    mu_crit = None
    for i in range(len(grid) - 1):
        if vals[i] >= 0.0 and vals[i + 1] <= 0.0:
            mu_crit = brentq(B0_corr, grid[i], grid[i + 1], xtol=1e-12)
            break
    if mu_crit is None:
        raise ValueError("No B_corr(0 deg) = 0 root found")

    target_abs_u_perp = Gamma_perp + math.log(1.0 + D_eff / delta)
    rhs = target_abs_u_perp / n_perp + 2.0 * _log1pexp(mu_crit - eps_1)
    eps_perp_int = 0.5 * (mu_crit - math.log(math.expm1(rhs)))
    return {
        "mu_l_crit": float(mu_crit),
        "eps_perp_int": float(eps_perp_int),
    }
