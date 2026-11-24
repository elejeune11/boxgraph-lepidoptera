import math
import numpy as np


def run(inputs, settings):
    R_0 = float(inputs["R_0"])
    u = float(inputs["u"])
    lam = float(inputs["lam"])
    gamma_I = float(inputs["gamma_I"])
    w = float(inputs["omega"])
    n_grid = int(settings.get("n_grid", 256))

    z = np.linspace(0.0, lam, n_grid, endpoint=False)
    dz = lam / float(n_grid)
    q = 2.0 * math.pi / lam
    delta_rho = u * np.sin(q * z) - (u**2) / (4.0 * R_0)
    rho = R_0 + delta_rho
    drho = u * q * np.cos(q * z)
    drho2 = drho**2
    sq = np.sqrt(1.0 + drho2)

    # Cancellation-free form of rho*sqrt(1 + drho^2) - R_0 = delta_rho + rho * drho^2 / (sq + 1)
    dA_um2 = float(2.0 * math.pi * np.sum(delta_rho + rho * drho2 / (sq + 1.0)) * dz)
    A0_um2 = 2.0 * math.pi * R_0 * lam
    A_um2 = A0_um2 + dA_um2
    Ao_um2 = float(2.0 * math.pi * np.sum(rho * drho2 / sq) * dz)

    # (nN/m) * (um^2) / um = (1e-9 N/m) * (1e-6 m) = 1e-15 N = 1e-3 pN
    F_over_lam_pN = gamma_I * (A_um2 + w * Ao_um2) / lam * 1e-3
    dF_over_lam_pN = gamma_I * (dA_um2 + w * Ao_um2) / lam * 1e-3
    return {
        "F_over_lam": F_over_lam_pN,
        "delta_F_over_lam": dF_over_lam_pN,
    }
