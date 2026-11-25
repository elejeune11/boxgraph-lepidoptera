import numpy as np
from scipy.optimize import minimize_scalar


def _F_tot(f_V, W_el, V_ratio):
    a1 = (max(0.0, f_V) * V_ratio) ** 0.2
    a2 = (max(0.0, 1.0 - f_V) * V_ratio) ** 0.2
    return (a1 ** 3) + (a2 ** 3) - W_el * (a1 + a2)


def run(inputs, settings):
    W_el = inputs["W_el"]
    V_ratio = inputs["V_ratio"]
    n_grid = max(10, int(round(float(settings.get("n_grid", 200)))))
    tol = float(settings.get("tol", 1e-8))

    # Log-spaced near 0 plus uniform grid to resolve sharp unequal minima near f_V -> 0
    n_log = max(5, n_grid // 4)
    n_lin = max(5, n_grid - n_log)
    grid = np.unique(np.concatenate([
        np.array([0.0]),
        np.geomspace(1e-8, 0.05, n_log),
        np.linspace(0.05, 0.5, n_lin),
    ]))
    vals = [_F_tot(float(f), W_el, V_ratio) for f in grid]
    idx = int(np.argmin(vals))
    lo = float(grid[max(0, idx - 1)])
    hi = float(grid[min(len(grid) - 1, idx + 1)])

    res = minimize_scalar(
        lambda f: _F_tot(f, W_el, V_ratio),
        bounds=(lo, hi),
        method="bounded",
        options={"xatol": tol},
    )
    candidates = [
        (0.0, _F_tot(0.0, W_el, V_ratio)),
        (0.5, _F_tot(0.5, W_el, V_ratio)),
        (float(res.x), float(res.fun)),
    ]
    best_f, best_F = min(candidates, key=lambda pair: pair[1])
    return {"f_V_eq": best_f, "F_min_dimless": best_F}
