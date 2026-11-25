import math
import numpy as np
from scipy.optimize import minimize_scalar


def _F_tot(f_V, W, V_ratio):
    x1 = (max(0.0, f_V) * V_ratio) ** (1.0 / 3.0)
    x2 = (max(0.0, 1.0 - f_V) * V_ratio) ** (1.0 / 3.0)
    f1 = (x1 ** 2) - W * (1.0 - 1.0 / math.sqrt(1.0 + x1 ** 2))
    f2 = (x2 ** 2) - W * (1.0 - 1.0 / math.sqrt(1.0 + x2 ** 2))
    return f1 + f2


def run(inputs, settings):
    W = inputs["W"]
    V_ratio = inputs["V_ratio"]
    n_grid = max(10, int(round(float(settings.get("n_grid", 200)))))
    tol = float(settings.get("tol", 1e-8))

    grid = np.linspace(0.0, 0.5, n_grid)
    vals = [_F_tot(float(f), W, V_ratio) for f in grid]
    idx = int(np.argmin(vals))
    lo = float(grid[max(0, idx - 1)])
    hi = float(grid[min(n_grid - 1, idx + 1)])

    res = minimize_scalar(
        lambda f: _F_tot(f, W, V_ratio),
        bounds=(lo, hi),
        method="bounded",
        options={"xatol": tol},
    )
    candidates = [
        (0.0, _F_tot(0.0, W, V_ratio)),
        (0.5, _F_tot(0.5, W, V_ratio)),
        (float(res.x), float(res.fun)),
    ]
    best_f, best_F = min(candidates, key=lambda pair: pair[1])
    return {"f_V_eq": best_f, "F_min_dimless": best_F}
