import numpy as np
from scipy.optimize import brentq


def run(inputs, settings):
    idx = int(round(float(inputs["idx"])))
    Dc = float(inputs["Dc"]); Da = float(inputs["Da"])
    A = float(inputs["A"]); Kd = float(inputs["Kd"])
    Om0 = float(inputs["Omega0"]); Om = float(inputs["Omega"])
    Omd0 = float(inputs["Omega_d0"]); Omd = float(inputs["Omega_d"])
    Z = float(inputs["Z"]); B = float(inputs["B"]); Nbar = float(inputs["Nbar"])

    ratio = A / Kd
    na_h = float(brentq(
        lambda na: Om0 * (1.0 + Om * na * na) * (Nbar - na) - (Omd0 + Omd * ratio * na) * na,
        0.0, Nbar, xtol=1e-14, rtol=1e-14,
    ))
    c_h = ratio * na_h
    ni_h = Nbar - na_h
    pi_p = 2.0 * Z * c_h - 3.0 * B * c_h * c_h

    J = np.zeros((8, 8), dtype=float)
    J[0, 4] = 1.0; J[1, 5] = 1.0; J[2, 6] = 1.0; J[3, 7] = 1.0
    J[4, 0] = Kd / Dc; J[4, 1] = -A / Dc; J[4, 7] = c_h / Dc
    dF_dC = Omd * na_h
    dF_dNa = -2.0 * Om0 * Om * na_h * ni_h + (Omd0 + Omd * c_h)
    dF_dNi = -Om0 * (1.0 + Om * na_h * na_h)
    J[5, 0] = dF_dC / Da; J[5, 1] = dF_dNa / Da; J[5, 2] = dF_dNi / Da; J[5, 7] = na_h / Da
    J[6, 0] = -dF_dC; J[6, 1] = -dF_dNa; J[6, 2] = -dF_dNi; J[6, 7] = ni_h
    J[7, 3] = 1.0; J[7, 4] = -pi_p

    evals = np.linalg.eigvals(J)
    # Collect the 4 distinct non-negative-real-part eigenvalues (one 0 + three Re > 0) ordered by (Re, Im)
    tol = 1e-6
    pos = [ev for ev in evals if ev.real > tol]
    pos.sort(key=lambda z: (round(float(z.real), 5), float(z.imag)))
    ordered = [0.0 + 0.0j] + pos  # 4 eigenvalues matching the 4 markers in Fig. 3(c)
    ev = ordered[max(0, min(len(ordered) - 1, idx - 1))]
    return {"re_part": float(ev.real), "im_part": float(ev.imag)}
