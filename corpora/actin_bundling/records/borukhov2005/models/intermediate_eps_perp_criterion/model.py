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

    # Printed Eq. 8 verbatim for gamma <= D/L: X = (|u_par| - Gamma_par)*L/D.
    # Note: because of the inverted sign, B0(mu_l) crosses zero upward from -inf to >0!
    def B0_printed(mu_l):
        u_par_tot = -n_par * _log1pexp(mu_l - 2.0 * eps_par) + 2.0 * n_par * _log1pexp(mu_l - eps_1)
        X = abs(u_par_tot) - Gamma_par * (L / D)
        if abs(X) < 1e-8:
            br = 0.5 * X
        elif X < -500.0:
            return -1e200
        else:
            br = 1.0 + (math.exp(-X) - 1.0) / X
        return t1 + c2 * br

    # Find crossing in the attractive regime mu_l < eps_1
    mu_hi = 2.0 * eps_par + 5.0
    mu_lo = 2.0 * eps_par - 20.0
    # Narrow bracket where B0_printed is finite
    grid = [mu_lo + (mu_hi - mu_lo) * i / 400.0 for i in range(401)]
    vals = [B0_printed(m) for m in grid]
    mu_crit = None
    for i in range(len(grid) - 1):
        if vals[i] * vals[i + 1] <= 0.0:
            mu_crit = brentq(B0_printed, grid[i], grid[i + 1], xtol=1e-12)
            break
    if mu_crit is None:
        raise ValueError("No B(0 deg) = 0 root found")

    # Solve B(90 deg) = 0 at mu_crit for eps_perp:
    # (8/pi) + (8*delta/(pi*D_eff))*(1 - exp(-Gamma_perp + |u_perp|)) = 0
    # => |u_perp| = Gamma_perp + ln(1 + D_eff/delta)
    target_abs_u_perp = Gamma_perp + math.log(1.0 + D_eff / delta)
    # In attractive regime u_perp < 0: |u_perp| = n_perp*ln(1+exp(mu-2*eps_perp)) - 2*n_perp*ln(1+exp(mu-eps_1))
    rhs = target_abs_u_perp / n_perp + 2.0 * _log1pexp(mu_crit - eps_1)
    eps_perp_int = 0.5 * (mu_crit - math.log(math.expm1(rhs)))
    return {
        "mu_l_crit": float(mu_crit),
        "eps_perp_int": float(eps_perp_int),
    }
