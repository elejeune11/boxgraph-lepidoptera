import math
import numpy as np
from scipy.optimize import brentq, minimize_scalar

_EULER_GAMMA = 0.5772156649015329
_C1_DISK = 128.0 / (45.0 * math.pi)
_C2_DISK = 0.01777119441631314


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
    h = kappa_inv / D_eff

    gs = np.linspace(D / L, 0.5 * math.pi, n_gamma)
    sg = np.sin(gs)
    base_B = (8.0 / math.pi) * sg * (1.0 - h * np.log(sg))

    c_cage = (math.pi / 4.0) * L * D_eff / (delta ** 2)
    sig_B = math.log(2.0 / ((D / L) ** 2))

    def mu_r_B(mu_l):
        u_par_tot = -n_par * _log1pexp(mu_l - 2.0 * eps_par) + 2.0 * n_par * _log1pexp(mu_l - eps_1)
        return math.log(c_cage) + sig_B + 3.0 * ((L / D) * Gamma_par + u_par_tot)

    def w_N_small(a):
        return (8.0 / math.pi) * (_C1_DISK * a * (1.0 + delta / D_eff - h * math.log(a)) - h * _C2_DISK * a)

    def mu_r_N(cr):
        def obj(la):
            a = math.exp(la)
            return cr * (math.log(2.0) - 2.0 * la) + 0.5 * (cr ** 2) * w_N_small(a)

        res = minimize_scalar(obj, bounds=(math.log(1e-4), math.log(1.0)), method="bounded")
        la = res.x
        a = math.exp(la)
        sig = math.log(2.0) - 2.0 * la
        w = w_N_small(a)
        return math.log(cr) + sig + cr * w

    mu_N_val = mu_r_N(c_r)

    def dilute_mu_and_cl(mu_l):
        up = -n_perp * _log1pexp(mu_l - 2.0 * eps_perp) + 2.0 * n_perp * _log1pexp(mu_l - eps_1)
        dup = n_perp * _sigmoid(mu_l - 2.0 * eps_perp) - 2.0 * n_perp * _sigmoid(mu_l - eps_1)
        B = base_B + (8.0 * delta / (math.pi * D_eff)) * sg * (1.0 - np.exp(-Gamma_perp / sg + abs(up)))
        w_I = float(np.trapezoid(B * sg, gs))
        mu_I = math.log(c_r) + c_r * w_I
        cl_free = (L ** 2 * D_eff / v_l) * math.exp(mu_l)
        cl_1rod = n_par * c_r * _sigmoid(mu_l - eps_1)
        if mu_I <= mu_N_val:
            dB = -(8.0 * delta / (math.pi * D_eff)) * sg * np.exp(-Gamma_perp / sg + abs(up)) * dup
            dw = float(np.trapezoid(dB * sg, gs))
            return mu_I, cl_free + cl_1rod - 0.5 * (c_r ** 2) * dw
        return mu_N_val, cl_free + cl_1rod

    def diff_fn(mu_l):
        return mu_r_B(mu_l) - dilute_mu_and_cl(mu_l)[0]

    # Scan upward from low mu_l to find the first bundling coexistence crossing
    grid = np.linspace(2.0 * eps_par - 10.0, 2.0 * eps_par + 2.0, 120)
    prev_m, prev_f = float(grid[0]), diff_fn(float(grid[0]))
    mu_eq = None
    for m in grid[1:]:
        mf = float(m)
        fv = diff_fn(mf)
        if prev_f * fv <= 0.0:
            mu_eq = brentq(diff_fn, prev_m, mf, xtol=1e-10)
            break
        prev_m, prev_f = mf, fv
    if mu_eq is None:
        raise ValueError("No bundle coexistence crossing found")
    _, c_l = dilute_mu_and_cl(mu_eq)
    return {"c_l": float(c_l), "mu_l_bundle": float(mu_eq)}
