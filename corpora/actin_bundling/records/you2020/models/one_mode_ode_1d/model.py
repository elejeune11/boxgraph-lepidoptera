import math
import numpy as np
from scipy.integrate import solve_ivp


def run(inputs, settings):
    aA = float(inputs["alpha_A"])
    aB = float(inputs["alpha_B"])
    k = float(inputs["kappa"])
    d = float(inputs["delta"])
    psi_probe = float(inputs["psi_probe"])

    rtol = float(settings.get("rtol", 1.0e-8))
    atol = float(settings.get("atol", 1.0e-11))
    t_max = float(settings.get("t_max", 25000.0))
    psi_init = float(settings.get("psi_init", 1.0e-4))

    # Start near F_S with a tiny positive phase perturbation theta = pi + psi_init
    arg_fs = (k * k - d * d - aA * aB) / aB
    if arg_fs > 0.0:
        rA0 = math.sqrt(arg_fs)
        rB0 = (k + d) * rA0 / aB
    else:
        rA0 = 0.05
        rB0 = 0.05
    th0 = math.pi + psi_init

    def rhs(t, y):
        rA, rB, th = y
        rA = max(rA, 1e-12)
        rB = max(rB, 1e-12)
        ct = math.cos(th)
        st = math.sin(th)
        drA = -(aA + rA * rA) * rA - (k - d) * rB * ct
        drB = -aB * rB - (k + d) * rA * ct
        dth = ((k - d) * rB / rA + (k + d) * rA / rB) * st
        return [drA, drB, dth]

    sol = solve_ivp(rhs, (0.0, t_max), [rA0, rB0, th0], method="RK45", rtol=rtol, atol=atol, dense_output=True)
    rA_f, rB_f, th_f = sol.y[:, -1]
    dPhi_f = ((k - d) * rB_f / rA_f - (k + d) * rA_f / rB_f) * math.sin(th_f)
    v_f = 0.5 * abs(dPhi_f)
    psi_f = abs(th_f - math.pi)
    amp_ratio_f = rA_f / rB_f
    v_ratio_f = -((k - d) * rB_f * rB_f) / ((k + d) * rA_f * rA_f)

    # Sample trajectory at psi = psi_probe along F_S -> F_T orbit
    psi_arr = sol.y[2] - math.pi
    v_ratio_probe = None
    amp_ratio_probe = None
    if psi_arr.max() >= psi_probe >= psi_arr.min():
        idx = int(np.searchsorted(psi_arr, psi_probe))
        if 0 < idx < len(psi_arr):
            t_lo, t_hi = sol.t[idx - 1], sol.t[idx]
            # Bisect on dense output for exact psi_probe
            for _ in range(40):
                t_mid = 0.5 * (t_lo + t_hi)
                ym = sol.sol(t_mid)
                if ym[2] - math.pi < psi_probe:
                    t_lo = t_mid
                else:
                    t_hi = t_mid
            yp = sol.sol(0.5 * (t_lo + t_hi))
            amp_ratio_probe = float(yp[0] / yp[1])
            v_ratio_probe = float(-((k - d) * yp[1] ** 2) / ((k + d) * yp[0] ** 2))

    return {
        "rho_A_final": float(rA_f),
        "rho_B_final": float(rB_f),
        "v_final": float(v_f),
        "psi_final": float(psi_f),
        "amp_ratio_final": float(amp_ratio_f),
        "v_ratio_final": float(v_ratio_f),
        "v_ratio_at_probe": v_ratio_probe,
        "amp_ratio_at_probe": amp_ratio_probe,
    }
