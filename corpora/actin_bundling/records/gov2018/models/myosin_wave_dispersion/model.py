import numpy as np


def run(inputs, settings):
    q = float(inputs["q"])
    kappa = float(inputs["kappa"])
    sigma = float(inputs["sigma"])
    kBT = float(inputs["kBT"])
    Lambda = float(inputs["Lambda"])
    H0 = float(inputs["H0"])
    ns = float(inputs["ns"])
    n0 = float(inputs["n0"])
    J = float(inputs["J"])
    xi = float(inputs["xi"])
    A = float(inputs["A"])
    A_star = float(inputs["A_star"])
    k_on = float(inputs["k_on"])
    k_off = float(inputs["k_off"])

    q2 = q * q
    q4 = q2 * q2
    D = Lambda * kBT
    Gh = (kappa * q4 + sigma * q2) / xi
    Gn = (
        D * q2
        + (kappa * Lambda * (H0 ** 2) * n0 / ns) * q2
        - (2.0 * J * Lambda * n0 / ns) * q2 * (1.0 - q2 / ns)
    )
    m12 = (A - kappa * H0 * q2) / xi
    m13 = -A_star / xi
    m21 = -(kappa * Lambda * H0 * n0 / ns) * q4

    M = np.array(
        [
            [-Gh, m12, m13],
            [m21, -Gn, 0.0],
            [0.0, k_on, -k_off],
        ],
        dtype=float,
    )
    eigs = np.linalg.eigvals(M)
    # Dominant mode: largest real part; for complex conjugate pairs take positive imaginary part
    dominant = max(eigs, key=lambda z: (round(float(z.real), 12), abs(float(z.imag))))
    w_real = float(dominant.real)
    w_imag = abs(float(dominant.imag))
    return {
        "omega_real": w_real,
        "omega_imag": w_imag,
        "phase_speed": w_imag / q if q > 0.0 else 0.0,
    }
