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

    # (nN/m) * (um^2) = 1e-21 J
    F_J = gamma_I * (A_um2 + w * Ao_um2) * 1e-21
    dF_J = gamma_I * (dA_um2 + w * Ao_um2) * 1e-21
    return {
        "A": A_um2,
        "A_omega": Ao_um2,
        "F": F_J,
        "delta_F": dF_J,
    }
