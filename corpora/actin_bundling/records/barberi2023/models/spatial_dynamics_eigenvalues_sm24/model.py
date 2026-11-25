import numpy as np
from scipy.optimize import brentq


def run(inputs, settings):
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
    Nbar = float(inputs["Nbar"])

    ratio = A / Kd
    na_h = float(
        brentq(
            lambda na: Om0 * (1.0 + Om * na * na) * (Nbar - na) - (Omd0 + Omd * ratio * na) * na,
            0.0,
            Nbar,
            xtol=1e-14,
            rtol=1e-14,
        )
    )
    c_h = ratio * na_h
    ni_h = Nbar - na_h
    pi_p = 2.0 * Z * c_h - 3.0 * B * c_h * c_h

    # Exact 8x8 Jacobian of Eqs. (SM24)-(SM31) at H for state [C, Na, Ni, V, PC, Pa, Pi, PV]
    J = np.zeros((8, 8), dtype=float)
    # Eqs. (SM24)-(SM27): dX(C, Na, Ni, V) = (PC, Pa, Pi, PV)
    J[0, 4] = 1.0
    J[1, 5] = 1.0
    J[2, 6] = 1.0
    J[3, 7] = 1.0
    # Eq. (SM28): dX PC = (1/Dc) * (V*PC + C*PV - A*Na + Kd*C)
    J[4, 0] = Kd / Dc
    J[4, 1] = -A / Dc
    J[4, 7] = c_h / Dc
    # Eq. (SM29): dX Pa = (1/Da) * [V*Pa + Na*PV - Om0*(1 + Om*Na^2)*Ni + (Omd0 + Omd*C)*Na]
    dF_dC = Omd * na_h
    dF_dNa = -2.0 * Om0 * Om * na_h * ni_h + (Omd0 + Omd * c_h)
    dF_dNi = -Om0 * (1.0 + Om * na_h * na_h)
    J[5, 0] = dF_dC / Da
    J[5, 1] = dF_dNa / Da
    J[5, 2] = dF_dNi / Da
    J[5, 7] = na_h / Da
    # Eq. (SM30): dX Pi = V*Pi + Ni*PV + Om0*(1 + Om*Na^2)*Ni - (Omd0 + Omd*C)*Na
    J[6, 0] = -dF_dC
    J[6, 1] = -dF_dNa
    J[6, 2] = -dF_dNi
    J[6, 7] = ni_h
    # Eq. (SM31): dX PV = V - (2*Z*C - 3*B*C^2)*PC
    J[7, 3] = 1.0
    J[7, 4] = -pi_p

    evals = np.linalg.eigvals(J)
    tol = 1e-6
    unit_mult = int(np.sum(np.abs(evals.real - 1.0) < tol))
    zero_mult = int(np.sum(np.abs(evals) < tol))
    max_re = float(np.max(evals.real))
    max_im = float(np.max(np.abs(evals.imag)))
    pos_re = evals.real[evals.real > tol]
    min_pos_re = float(np.min(pos_re)) if len(pos_re) else 0.0
    has_complex = 1.0 if max_im > tol else 0.0

    return {
        "max_real_part": max_re,
        "max_imag_part": max_im,
        "min_nonzero_real_part": min_pos_re,
        "unit_eval_multiplicity": float(unit_mult),
        "zero_eval_multiplicity": float(zero_mult),
        "has_complex_eigenvalues": has_complex,
    }
