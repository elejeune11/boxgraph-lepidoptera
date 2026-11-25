import math
import numpy as np
from scipy.optimize import brentq


def run(inputs, settings):
    N_loc = float(inputs["N_loc"])
    W = float(inputs["W"])
    Dc = float(inputs["Dc"])
    Da = float(inputs["Da"])
    A = float(inputs["A"])
    Kd = float(inputs["Kd"])
    Om0 = float(inputs["Omega0"])
    Om = float(inputs["Omega"])
    Omd0 = float(inputs["Omega_d0"])
    Omd = float(inputs["Omega_d"])
    Z = float(inputs["Z"])
    B = float(inputs["B"])

    n_modes = int(round(settings.get("n_modes", 120)))
    dK = math.pi / W
    Ks = np.arange(1, n_modes + 1, dtype=float) * dK

    ratio = A / Kd
    na_h = float(
        brentq(
            lambda na: Om0 * (1.0 + Om * na * na) * (N_loc - na) - (Omd0 + Omd * ratio * na) * na,
            0.0,
            N_loc,
            xtol=1e-14,
            rtol=1e-14,
        )
    )
    c_h = ratio * na_h
    nh_m = 2.0 * na_h - N_loc
    pi_p = 2.0 * Z * c_h - 3.0 * B * c_h * c_h

    k2 = Ks * Ks
    fac = (k2 / (1.0 + k2)) * pi_p
    da_p = Da + 1.0
    da_m = Da - 1.0
    n = len(Ks)
    # Evaluates printed Eqs. (SM15)-(SM23) verbatim (with N_H^+ = 1 left unscaled)
    M = np.empty((n, 3, 3), dtype=float)
    M[:, 0, 0] = fac * c_h - Dc * k2 - Kd
    M[:, 0, 1] = 0.5 * A
    M[:, 0, 2] = 0.5 * A
    M[:, 1, 0] = fac
    M[:, 1, 1] = -0.5 * da_p * k2
    M[:, 1, 2] = -0.5 * da_m * k2
    M[:, 2, 0] = fac * nh_m - (1.0 + nh_m) * Omd
    M[:, 2, 1] = (
        -0.5 * da_m * k2
        + Om0
        - 0.25 * (nh_m - 3.0) * (nh_m + 1.0) * Om * Om0
        - Omd0
        - c_h * Omd
    )
    M[:, 2, 2] = (
        -0.5 * da_p * k2
        - 0.25 * (4.0 + (1.0 + nh_m) * (3.0 * nh_m - 1.0) * Om) * Om0
        - Omd0
        - c_h * Omd
    )

    evals = np.linalg.eigvals(M)
    growth = np.max(evals.real, axis=1)
    best_i = int(np.argmax(growth))
    unst = Ks[growth > 0.0]
    if len(unst) == 0:
        k_left = float("nan")
        k_right = float("nan")
    else:
        k_left = float(unst[0])
        k_right = float(unst[-1])

    return {
        "K_star": float(Ks[best_i]),
        "K_left": k_left,
        "K_right": k_right,
        "delta_K": float(dK),
        "max_growth_rate": float(growth[best_i]),
    }
