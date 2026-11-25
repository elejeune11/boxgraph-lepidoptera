import numpy as np
from scipy.optimize import brentq


def run(inputs, settings):
    K = float(inputs["K"])
    Dc = float(inputs["Dc"])
    Da = float(inputs["Da"])
    A = float(inputs["A"])
    Kd = float(inputs["Kd"])
    Om0 = float(inputs["Omega0"])
    Om = float(inputs["Omega"])
    Omd0 = float(inputs["Omega_d0"])
    Omd = float(inputs["Omega_d"])

    ratio = A / Kd
    na_h = float(
        brentq(
            lambda na: Om0 * (1.0 + Om * na * na) * (1.0 - na) - (Omd0 + Omd * ratio * na) * na,
            0.0,
            1.0,
            xtol=1e-14,
            rtol=1e-14,
        )
    )
    c_h = ratio * na_h
    nh_m = 2.0 * na_h - 1.0

    k2 = K * K
    da_p = Da + 1.0
    da_m = Da - 1.0

    M = np.array(
        [
            [-Dc * k2 - Kd, 0.5 * A, 0.5 * A],
            [0.0, -0.5 * da_p * k2, -0.5 * da_m * k2],
            [
                -(1.0 + nh_m) * Omd,
                -0.5 * da_m * k2
                + Om0
                - 0.25 * (nh_m - 3.0) * (nh_m + 1.0) * Om * Om0
                - Omd0
                - c_h * Omd,
                -0.5 * da_p * k2
                - 0.25 * (4.0 + (1.0 + nh_m) * (3.0 * nh_m - 1.0) * Om) * Om0
                - Omd0
                - c_h * Omd,
            ],
        ],
        dtype=float,
    )
    evals = np.linalg.eigvals(M)
    idx = int(np.argmax(evals.real))
    lead = evals[idx]
    return {
        "max_growth_rate": float(lead.real),
        "dominant_freq": float(abs(lead.imag)),
    }
