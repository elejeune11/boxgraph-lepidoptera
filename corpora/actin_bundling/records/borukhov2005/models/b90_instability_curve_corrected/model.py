import math
import numpy as np
from scipy.optimize import brentq

_EULER_GAMMA = 0.5772156649015329


def _log1pexp(x):
    if x > 50.0:
        return x
    if x < -50.0:
        return math.exp(x)
    return math.log1p(math.exp(x))


def _sigmoid(x):
    if x >= 0.0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


def run(inputs, settings):
    c_r = inputs["c_r"]
    D = inputs["D"]
    L = inputs["L"]
    delta = inputs["delta"]
    kappa_inv = inputs["kappa_inv"]
    n_perp = inputs["n_perp"]
    Gamma_par = inputs["Gamma_par"]
    Gamma_perp = inputs["Gamma_perp"]
    eps_1 = inputs["eps_1"]
    eps_par = inputs["eps_par"]
    eps_perp = inputs["eps_perp"]
    v_l = inputs["v_l"]
    n_gamma = int(settings.get("n_gamma", 2000))

    D_eff = D + kappa_inv * (_EULER_GAMMA + math.log(Gamma_perp))
    n_par = L / delta

    target_abs_u_perp = Gamma_perp + math.log(1.0 + D_eff / delta)

    def u_perp_fn(mu):
        return -n_perp * _log1pexp(mu - 2.0 * eps_perp) + 2.0 * n_perp * _log1pexp(mu - eps_1)

    mu_star = brentq(lambda m: -u_perp_fn(m) - target_abs_u_perp, 2.0 * eps_perp - 10.0, eps_1 + 2.0, xtol=1e-12)
    up_star = u_perp_fn(mu_star)
    d_abs_up_dmu = n_perp * _sigmoid(mu_star - 2.0 * eps_perp) - 2.0 * n_perp * _sigmoid(mu_star - eps_1)

    cl_free = (L ** 2 * D_eff / v_l) * math.exp(mu_star)
    cl_1rod = n_par * c_r * _sigmoid(mu_star - eps_1)

    gammas = np.linspace(D / L, 0.5 * math.pi, n_gamma)
    sg = np.sin(gammas)
    dB_large_dmu = -(8.0 * delta / (math.pi * D_eff)) * sg * np.exp(-Gamma_perp / sg + abs(up_star)) * d_abs_up_dmu
    dw_large_dmu = float(np.trapezoid(dB_large_dmu * sg, gammas))

    # Corrected Eq. 8 for gamma <= D/L: Y = (Gamma_par - |u_par|) * L/D
    u_par_tot = -n_par * _log1pexp(mu_star - 2.0 * eps_par) + 2.0 * n_par * _log1pexp(mu_star - eps_1)
    d_abs_upar_tot_dmu = n_par * _sigmoid(mu_star - 2.0 * eps_par) - 2.0 * n_par * _sigmoid(mu_star - eps_1)
    Y = Gamma_par * (L / D) - abs(u_par_tot)
    if abs(Y) < 1e-6:
        dbr_dY = 0.5
    else:
        dbr_dY = (1.0 - (1.0 + Y) * math.exp(-Y)) / (Y ** 2)
    dB_small_dmu = -16.0 * (D * delta / (D_eff * L)) * dbr_dY * d_abs_upar_tot_dmu
    cap_weight = 1.0 - math.cos(D / L)
    dw_small_dmu = cap_weight * dB_small_dmu

    c_l = cl_free + cl_1rod - 0.5 * (c_r ** 2) * (dw_small_dmu + dw_large_dmu)
    return {"c_l": float(c_l), "mu_l_star": float(mu_star)}
