import cmath
import math
import numpy as np
from scipy.optimize import minimize_scalar

def run(inputs, settings):
    kd = inputs["k_d"]
    L = inputs["L"]
    c0 = inputs["c_0"]
    b = inputs["b"]
    er = inputs["eta_rot"]
    beta = inputs["beta"]
    lo = inputs["lambda_o"]
    kap = inputs["kappa"]
    q0 = inputs["q_0"]
    ratio = inputs["activity_ratio"]

    G = c0 if c0 >= 0.0 else 8.0 * b * q0 * q0

    def lam_fn(x):
        if x <= 0.0:
            return 1e30
        num = (1.0 / x + 4.0) * (x + kd) * (2.0 * G + 2.0 * L * x)
        den = (1.0 + kap * q0) * (2.0 * G + 2.0 * L * x) + kap * (beta * (x + kd) + 2.0 * lo * q0)
        return num / den if den > 0.0 else 1e30

    n_grid = int(settings.get("n_grid", 2000))
    xs = np.geomspace(1e-5, 1e4, n_grid)
    vals = np.array([lam_fn(float(x)) for x in xs])
    i0 = int(np.argmin(vals))
    res_c = minimize_scalar(
        lam_fn,
        bounds=(float(xs[max(0, i0 - 3)]), float(xs[min(n_grid - 1, i0 + 3)])),
        method="bounded",
        options={"xatol": 1e-13},
    )
    lam_c = float(res_c.fun)
    nu_c = math.sqrt(float(res_c.x))
    lam = ratio * lam_c

    def alpha_plus(nu):
        x = nu * nu
        f_vis = 4.0 + 1.0 / x
        g_nem = 2.0 * G + 2.0 * L * x
        A = 2.0 * er * f_vis - beta * beta
        B = (
            f_vis * g_nem
            + 2.0 * er * (f_vis * (x + kd) - lam * (1.0 + kap * q0))
            - beta * (beta * (x + kd) + 2.0 * lo * q0 + lam * kap)
        )
        C = (f_vis * (x + kd) - lam * (1.0 + kap * q0)) * g_nem - lam * kap * (
            beta * (x + kd) + 2.0 * lo * q0
        )
        disc = cmath.sqrt(B * B - 4.0 * A * C)
        r1 = ((-B + disc) / (2.0 * A)).real
        r2 = ((-B - disc) / (2.0 * A)).real
        return max(r1, r2)

    nus = np.geomspace(1e-3, 1e2, n_grid)
    gvals = np.array([alpha_plus(float(nu)) for nu in nus])
    im = int(np.argmax(gvals))
    res_m = minimize_scalar(
        lambda nu: -alpha_plus(nu),
        bounds=(float(nus[max(0, im - 3)]), float(nus[min(n_grid - 1, im + 3)])),
        method="bounded",
        options={"xatol": 1e-12},
    )
    nu_max = float(res_m.x)
    alpha_max = float(-res_m.fun)
    return {
        "lambda_crit": lam_c,
        "nu_crit": nu_c,
        "nu_max": nu_max,
        "alpha_max": alpha_max,
        "ell_max": 2.0 * math.pi / nu_max,
        "wavelength_ratio": nu_c / nu_max,
    }
