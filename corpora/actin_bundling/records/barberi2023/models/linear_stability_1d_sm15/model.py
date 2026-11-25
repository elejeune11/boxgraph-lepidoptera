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
    Z = float(inputs["Z"])
    B = float(inputs["B"])

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
    pi_p = 2.0 * Z * c_h - 3.0 * B * c_h * c_h

    k2 = K * K
    fac = k2 / (1.0 + k2) * pi_p
    da_p = Da + 1.0
    da_m = Da - 1.0

    # Exact printed Eqs. (SM15)-(SM23) for N_H^+ = 1
    M = np.array(
        [
            [fac * c_h - Dc * k2 - Kd, 0.5 * A, 0.5 * A],
            [fac, -0.5 * da_p * k2, -0.5 * da_m * k2],
            [
                fac * nh_m - (1.0 + nh_m) * Omd,
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
        "vel_response_imag": float((K / (1.0 + k2)) * pi_p),
    }
