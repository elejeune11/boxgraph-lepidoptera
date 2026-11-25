import numpy as np
from scipy.optimize import brentq


def run(inputs, settings):
    Z = float(inputs["Z"])
    B = Z  # SM Table II sets B = Z
    Dc = float(inputs["Dc"])
    Da = float(inputs["Da"])
    A = float(inputs["A"])
    Kd = float(inputs["Kd"])
    Om0 = float(inputs["Omega0"])
    Omd0 = float(inputs["Omega_d0"])
    Omd = float(inputs["Omega_d"])

    n_k = int(round(settings.get("n_k", 400)))
    k_min = float(settings.get("k_min", 0.02))
    k_max = float(settings.get("k_max", 10.0))
    omega_max = float(settings.get("omega_max", 40.0))

    # Log-linear combined grid resolves both long-wave K -> 0 and finite-K peaks
    n_low = max(20, n_k // 4)
    n_hi = max(20, n_k - n_low)
    Ks = np.unique(
        np.concatenate(
            [
                np.geomspace(k_min, 0.5, n_low),
                np.linspace(0.5, k_max, n_hi),
            ]
        )
    )
    k2 = Ks * Ks
    k_fac = k2 / (1.0 + k2)
    da_p = Da + 1.0
    da_m = Da - 1.0
    ratio = A / Kd

    def max_growth(Om):
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

        # Vectorized 3x3 eigenvalues over all K
        n = len(Ks)
        M = np.empty((n, 3, 3), dtype=float)
        M[:, 0, 0] = k_fac * pi_p * c_h - Dc * k2 - Kd
        M[:, 0, 1] = 0.5 * A
        M[:, 0, 2] = 0.5 * A
        M[:, 1, 0] = k_fac * pi_p
        M[:, 1, 1] = -0.5 * da_p * k2
        M[:, 1, 2] = -0.5 * da_m * k2
        M[:, 2, 0] = k_fac * pi_p * nh_m - (1.0 + nh_m) * Omd
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
        # Leading eigenvalue per K
        best_col = np.argmax(evals.real, axis=1)
        lead_per_k = evals[np.arange(n), best_col]
        best_i = int(np.argmax(lead_per_k.real))
        return float(lead_per_k[best_i].real), float(Ks[best_i]), float(abs(lead_per_k[best_i].imag))

    g0, k0, f0 = max_growth(0.0)
    if g0 >= 0.0:
        om_c, k_c, f_c = 0.0, k0, f0
    else:
        om_c = float(brentq(lambda om: max_growth(om)[0], 0.0, omega_max, xtol=1e-9, rtol=1e-9))
        _, k_c, f_c = max_growth(om_c + 1e-5)

    return {
        "Omega_crit": om_c,
        "K_crit": k_c,
        "freq_crit": f_c,
        "is_oscillatory": 1.0 if f_c > 1e-4 else 0.0,
    }
