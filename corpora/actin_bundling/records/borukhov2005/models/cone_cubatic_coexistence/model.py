import math
import numpy as np
from scipy.optimize import brentq, minimize_scalar

_EULER_GAMMA = 0.5772156649015329
_C1_DISK = 128.0 / (45.0 * math.pi)
# Exact small-angle disk average <(2u) ln(2u)> where P(u) = (16/pi)*u*(arccos(u) - u*sqrt(1-u^2))
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
    c_l_target = inputs["c_l"]
    D = inputs["D"]
    L = inputs["L"]
    delta = inputs["delta"]
    kappa_inv = inputs["kappa_inv"]
    n_perp = inputs["n_perp"]
    Gamma_perp = inputs["Gamma_perp"]
    eps_1 = inputs["eps_1"]
    eps_perp = inputs["eps_perp"]
    v_l = inputs["v_l"]
    n_gamma = int(settings.get("n_gamma", 2000))

    D_eff = D + kappa_inv * (_EULER_GAMMA + math.log(Gamma_perp))
    n_par = L / delta
    h = kappa_inv / D_eff

    gs = np.linspace(D / L, 0.5 * math.pi, n_gamma)
    sg = np.sin(gs)
    base_B = (8.0 / math.pi) * sg * (1.0 - h * np.log(sg))

    def u_perp_and_deriv(mu_l):
        up = -n_perp * _log1pexp(mu_l - 2.0 * eps_perp) + 2.0 * n_perp * _log1pexp(mu_l - eps_1)
        d_abs_up = n_perp * _sigmoid(mu_l - 2.0 * eps_perp) - 2.0 * n_perp * _sigmoid(mu_l - eps_1)
        return up, d_abs_up

    def B90(up):
        return (8.0 / math.pi) * (1.0 + (delta / D_eff) * (1.0 - math.exp(-Gamma_perp + abs(up))))

    def w_I(up):
        B = base_B + (8.0 * delta / (math.pi * D_eff)) * sg * (1.0 - np.exp(-Gamma_perp / sg + abs(up)))
        return float(np.trapezoid(B * sg, gs))

    def w_N_small(a):
        return (8.0 / math.pi) * (_C1_DISK * a * (1.0 + delta / D_eff - h * math.log(a)) - h * _C2_DISK * a)

    def opt_N(cr):
        def obj(la):
            a = math.exp(la)
            return cr * (math.log(2.0) - 2.0 * la) + 0.5 * (cr ** 2) * w_N_small(a)

        res = minimize_scalar(obj, bounds=(math.log(1e-4), math.log(1.0)), method="bounded")
        la = res.x
        a = math.exp(la)
        sig = math.log(2.0) - 2.0 * la
        w = w_N_small(a)
        return math.log(cr) + sig + cr * w, cr + 0.5 * (cr ** 2) * w

    def opt_C(cr, b90):
        def obj(la):
            a = math.exp(la)
            w = (1.0 / 3.0) * w_N_small(a) + (2.0 / 3.0) * b90
            return cr * (math.log(2.0 / 3.0) - 2.0 * la) + 0.5 * (cr ** 2) * w

        res = minimize_scalar(obj, bounds=(math.log(1e-4), math.log(0.25 * math.pi)), method="bounded")
        la = res.x
        a = math.exp(la)
        sig = math.log(2.0 / 3.0) - 2.0 * la
        w = (1.0 / 3.0) * w_N_small(a) + (2.0 / 3.0) * b90
        return math.log(cr) + sig + cr * w, cr + 0.5 * (cr ** 2) * w

    def cl_in_C(cC, mu_l):
        up, d_abs_up = u_perp_and_deriv(mu_l)
        cl_free = (L ** 2 * D_eff / v_l) * math.exp(mu_l)
        cl_1rod = n_par * cC * _sigmoid(mu_l - eps_1)
        cl_2body = 0.5 * (cC ** 2) * (2.0 / 3.0) * (8.0 * delta / (math.pi * D_eff)) * math.exp(-Gamma_perp + abs(up)) * d_abs_up
        return cl_free + cl_1rod + cl_2body

    # Upper limit of mu_l where B(90 deg) = 0
    target_abs_up = Gamma_perp + math.log(1.0 + D_eff / delta)
    mu_b90 = brentq(lambda m: -u_perp_and_deriv(m)[0] - target_abs_up, 2.0 * eps_perp - 10.0, eps_1 + 2.0, xtol=1e-12)

    # 1) C/(N+C) right boundary: find mu_l where C-N coexistence has cl_C(mu_l) = c_l_target
    def cC_CN_at_mu(mu_l):
        up, _ = u_perp_and_deriv(mu_l)
        b90 = B90(up)

        def diff_CN(cC):
            muC, PC = opt_C(cC, b90)
            hi_n = max(100.0, PC * 1.5)
            while opt_N(hi_n)[1] < PC:
                hi_n *= 2.0
            cN = brentq(lambda c: opt_N(c)[1] - PC, 0.1, hi_n, xtol=1e-7)
            muN, _ = opt_N(cN)
            return muC - muN

        # Scan for sign change of diff_CN(cC) on [3, 60]
        grid_c = [3.0 + i * 1.5 for i in range(40)]
        prev_c, prev_v = grid_c[0], diff_CN(grid_c[0])
        for c in grid_c[1:]:
            v = diff_CN(c)
            if prev_v * v <= 0.0:
                return brentq(diff_CN, prev_c, c, xtol=1e-6)
            prev_c, prev_v = c, v
        return None

    # Bracket mu_l for C-N where cl_in_C == c_l_target
    mu_grid = np.linspace(mu_b90 - 0.85, mu_b90 - 0.10, 25)
    prev_m, prev_f = None, None
    mu_cn = None
    for m in mu_grid:
        c_val = cC_CN_at_mu(float(m))
        if c_val is None:
            continue
        f_val = cl_in_C(c_val, float(m)) - c_l_target
        if prev_f is not None and prev_f * f_val <= 0.0:
            mu_cn = brentq(lambda mm: cl_in_C(cC_CN_at_mu(mm), mm) - c_l_target, prev_m, float(m), xtol=1e-7)
            break
        prev_m, prev_f = float(m), f_val
    c_C_right = cC_CN_at_mu(mu_cn)

    # 2) (I+C)/C left boundary: find mu_l where I-C coexistence has cl_C(mu_l) = c_l_target
    def cC_IC_at_mu(mu_l):
        up, _ = u_perp_and_deriv(mu_l)
        b90 = B90(up)
        wIv = w_I(up)

        def diff_IC(cC):
            muC, PC = opt_C(cC, b90)
            cI = (-1.0 + math.sqrt(1.0 + 2.0 * wIv * PC)) / wIv
            return muC - (math.log(cI) + cI * wIv)

        return brentq(diff_IC, 5.0, 60.0, xtol=1e-6)

    mu_ic = brentq(lambda m: cl_in_C(cC_IC_at_mu(m), m) - c_l_target, 2.0 * eps_perp - 4.0, mu_b90 - 0.1, xtol=1e-7)
    c_C_left = cC_IC_at_mu(mu_ic)

    return {"c_C_left": float(c_C_left), "c_C_right": float(c_C_right)}
