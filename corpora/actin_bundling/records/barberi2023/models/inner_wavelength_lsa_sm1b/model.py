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

    # Exact linearization of Eqs. (SM1)-(SM5) around H with total nucleator density N_bar = N_loc
    # (mathematically equivalent to rescaling Omega -> Omega*N_loc^2, Omega_d -> Omega_d*N_loc,
    # Z -> Z*N_loc^2, B -> B*N_loc^3 per SM Table II while keeping N_H^+ = 1).
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
    ni_h = N_loc - na_h
    pi_p = 2.0 * Z * c_h - 3.0 * B * c_h * c_h

    dRa_dC = -Omd * na_h
    dRa_dNa = 2.0 * Om0 * Om * na_h * ni_h - (Omd0 + Omd * c_h)
    dRa_dNi = Om0 * (1.0 + Om * na_h * na_h)

    k2 = Ks * Ks
    fac = (k2 / (1.0 + k2)) * pi_p
    n = len(Ks)
    M = np.empty((n, 3, 3), dtype=float)
    M[:, 0, 0] = fac * c_h - Dc * k2 - Kd
    M[:, 0, 1] = A
    M[:, 0, 2] = 0.0
    M[:, 1, 0] = fac * na_h + dRa_dC
    M[:, 1, 1] = -Da * k2 + dRa_dNa
    M[:, 1, 2] = dRa_dNi
    M[:, 2, 0] = fac * ni_h - dRa_dC
    M[:, 2, 1] = -dRa_dNa
    M[:, 2, 2] = -k2 - dRa_dNi

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
