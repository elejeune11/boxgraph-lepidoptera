import math
import numpy as np
from scipy.optimize import minimize_scalar

def run(inputs, settings):
    kd = inputs["k_d"]
    L = inputs["L"]
    c0 = inputs["c_0"]
    b = inputs["b"]
    beta = inputs["beta"]
    lo = inputs["lambda_o"]
    kap = inputs["kappa"]
    q0 = inputs["q_0"]

    G = c0 if c0 >= 0.0 else 8.0 * b * q0 * q0

    def lam_fn(x):
        if x <= 0.0:
            return 1e30
        num = (1.0 / x + 4.0) * (x + kd) * (2.0 * G + 2.0 * L * x)
        den = (1.0 + kap * q0) * (2.0 * G + 2.0 * L * x) + kap * (beta * (x + kd) + 2.0 * lo * q0)
        if den <= 0.0:
            return 1e30
        return num / den

    n_grid = int(settings.get("n_grid", 4000))
    x_min = float(settings.get("x_min", 1e-5))
    x_max = float(settings.get("x_max", 1e4))
    xs = np.geomspace(x_min, x_max, n_grid)
    vals = np.array([lam_fn(float(x)) for x in xs])
    idx = int(np.argmin(vals))
    lo_b = float(xs[max(0, idx - 3)])
    hi_b = float(xs[min(n_grid - 1, idx + 3)])
    res = minimize_scalar(lam_fn, bounds=(lo_b, hi_b), method="bounded", options={"xatol": 1e-13})
    x_c = float(res.x)
    lam_c = float(res.fun)
    nu_c = math.sqrt(x_c)
    return {
        "lambda_crit": lam_c,
        "nu_crit_sq": x_c,
        "nu_crit": nu_c,
        "ell_pattern": 2.0 * math.pi / nu_c,
    }
