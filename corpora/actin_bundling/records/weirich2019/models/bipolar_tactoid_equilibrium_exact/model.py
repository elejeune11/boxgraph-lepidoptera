import math
import numpy as np
from scipy.optimize import minimize_scalar


def _free_energy(alpha, K, gamma_pN, V):
    sa = math.sin(alpha)
    ca = math.cos(alpha)
    f_el = -math.pi * (7.0 * alpha * ca - 7.0 * sa + (alpha ** 2) * sa)
    f_s = 4.0 * math.pi * (sa - alpha * ca)
    v = 2.0 * math.pi * (sa - alpha * ca - (sa ** 3) / 3.0)
    R = (V / v) ** (1.0 / 3.0)
    return K * R * f_el + gamma_pN * (R ** 2) * f_s, R


def run(inputs, settings):
    K = inputs["K"]
    gamma_pN = inputs["gamma"] * 1e-3  # nN/m -> pN/um
    V = inputs["V"]
    n_grid = int(settings.get("n_grid", 200))
    tol = float(settings.get("tol", 1e-8))

    grid = np.linspace(1e-3, 0.5 * math.pi, n_grid)
    vals = [_free_energy(float(a), K, gamma_pN, V)[0] for a in grid]
    idx = int(np.argmin(vals))
    lo = float(grid[max(0, idx - 1)])
    hi = float(grid[min(n_grid - 1, idx + 1)])

    res = minimize_scalar(
        lambda a: _free_energy(a, K, gamma_pN, V)[0],
        bounds=(lo, hi),
        method="bounded",
        options={"xatol": tol},
    )
    alpha_eq = float(res.x)
    # Check right endpoint alpha = pi/2 explicitly (spherical saturation limit)
    F_eq, R_eq = _free_energy(alpha_eq, K, gamma_pN, V)
    F_sph, R_sph = _free_energy(0.5 * math.pi, K, gamma_pN, V)
    if F_sph <= F_eq:
        alpha_eq = 0.5 * math.pi
        F_eq, R_eq = F_sph, R_sph

    AR_eq = 1.0 / math.tan(0.5 * alpha_eq)
    return {
        "alpha_eq": alpha_eq,
        "R_eq": R_eq,
        "AR_eq": AR_eq,
        "F_eq": F_eq,
    }
