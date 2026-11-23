import math
import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import brentq


def run(inputs, settings):
    g = inputs["g"]
    p = inputs["p"]
    Y = inputs["Y"]
    nu = inputs["nu"]
    h = inputs["h"]
    L = inputs["L"]
    n_grid = int(round(settings.get("n_grid", 64)))
    tol = float(settings.get("tol", 1.0e-9))

    if p == 0.0:
        g_c = (math.pi ** 2) / (3.0 * (1.0 - nu ** 2)) * ((h / L) ** 2)
        if g <= g_c:
            return {"w_max": 0.0, "w_min": 0.0, "c": 0.0}
        a = (2.0 * h / math.sqrt(3.0)) * math.sqrt(g / g_c - 1.0)
        return {"w_max": a, "w_min": 0.0, "c": Y * (g - g_c)}

    D = Y * (h ** 3) / (12.0 * (1.0 - nu ** 2))
    Q_dim = p * ((L / 2.0) ** 4) / D
    fac = 12.0 * (1.0 - nu ** 2) * ((L / (2.0 * h)) ** 2)

    dxi = 2.0 / n_grid
    M = n_grid - 1
    rhs = np.full(M, Q_dim, dtype=float)
    inv_dxi4 = 1.0 / (dxi ** 4)
    inv_dxi2 = 1.0 / (dxi ** 2)
    # Exact discrete first even buckling eigenvalue for 2nd-order centered FD with ghost-point clamped BCs
    alpha_crit = (2.0 * math.sin(math.pi * dxi / 2.0) / dxi) ** 2

    def solve_for_alpha(alpha):
        ab = np.zeros((5, M), dtype=float)
        ab[0, 2:] = inv_dxi4
        ab[1, 1:] = -4.0 * inv_dxi4 + alpha * inv_dxi2
        ab[2, :] = 6.0 * inv_dxi4 - 2.0 * alpha * inv_dxi2
        ab[2, 0] += inv_dxi4
        ab[2, -1] += inv_dxi4
        ab[3, :-1] = -4.0 * inv_dxi4 + alpha * inv_dxi2
        ab[4, :-2] = inv_dxi4
        W_int = solve_banded((2, 2), ab, rhs)
        W = np.zeros(n_grid + 1, dtype=float)
        W[1:-1] = W_int
        dW = np.diff(W)
        I = float(np.sum(dW ** 2) / dxi)
        return W, I

    def res_log(log_delta):
        alpha = alpha_crit - math.exp(log_delta)
        _, I = solve_for_alpha(alpha)
        c_over_Y = I / ((1.0 - nu ** 2) * (L ** 2))
        return alpha - fac * (g - c_over_Y)

    ld_lo = -20.0
    while res_log(ld_lo) <= 0.0 and ld_lo > -60.0:
        ld_lo -= 5.0
    ld_hi = 4.0
    while res_log(ld_hi) >= 0.0 and ld_hi < 40.0:
        ld_hi += 5.0

    ld_sol = brentq(res_log, ld_lo, ld_hi, rtol=tol, xtol=1.0e-12, maxiter=500)
    alpha_sol = alpha_crit - math.exp(ld_sol)
    W, I = solve_for_alpha(alpha_sol)
    c_val = Y * I / ((1.0 - nu ** 2) * (L ** 2))
    return {
        "w_max": float(np.max(W)),
        "w_min": float(np.min(W)),
        "c": float(c_val),
    }
