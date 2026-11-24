import math
import numpy as np
from scipy.optimize import brentq


def run(inputs, settings):
    phibar = inputs["phibar"]
    k_n = inputs["k_norm"]
    kb_n = inputs["k_b_norm"]
    psi = inputs["psi"]
    Vs_n = inputs["Vs_norm"]
    Q_n = inputs["Q_norm"]
    w = inputs["w"]
    D_B = inputs["D_B"]

    k_0 = D_B / (w ** 2)
    k = k_n * k_0
    k_b = kb_n * k_0
    k_s = k * psi
    gb = w / 6.0
    V_s = Vs_n * (w ** 3)
    Q = Q_n * (w ** 3) * k_0

    R_crit_inf = None
    R_bar_eq23 = None
    R_bar_eqB1 = None
    if phibar > k_b / k + 1e-12:
        R_crit_inf = (gb / (phibar - k_b / k)) / w
        V0 = (V_s / psi) * (phibar - k_b / k)
        R_bar_eq23 = ((3.0 * V0 / (4.0 * math.pi)) ** (1.0 / 3.0)) / w
        dV_b1 = (Q * (k_b ** 2) - 4.0 * math.pi * D_B * (k ** 2) * phibar * gb) / ((k * phibar - k_b) * k * k_b * psi)
        if V0 + dV_b1 > 0.0:
            R_bar_eqB1 = ((3.0 * (V0 + dV_b1) / (4.0 * math.pi)) ** (1.0 / 3.0)) / w

    A = k * (psi ** 2)
    B = -V_s * psi * (k * phibar - k_b - Q / V_s)
    C = -V_s * (Q * phibar - 4.0 * math.pi * D_B * gb)
    disc = B ** 2 - 4.0 * A * C
    R_bar_quad = None
    if disc >= 0.0:
        Vq = (-B + math.sqrt(disc)) / (2.0 * A)
        if Vq > 0.0:
            R_bar_quad = ((3.0 * Vq / (4.0 * math.pi)) ** (1.0 / 3.0)) / w

    # Full Eq. (10) + Eq. (12) single droplet
    def dV_full(R):
        V = (4.0 * math.pi / 3.0) * (R ** 3)
        if V >= 0.999999 * V_s:
            return None
        A0 = phibar - (V / V_s) * psi
        cB = (V_s - V) / V_s
        qa = -k * cB
        qb = -k_b + k * A0 - 4.0 * math.pi * D_B * R / (V_s - V)
        qc = 4.0 * math.pi * D_B * gb / (V_s - V)
        dc = qb ** 2 - 4.0 * qa * qc
        if dc < 0.0:
            return None
        pB = max(0.0, (-qb - math.sqrt(dc)) / (2.0 * qa))
        pA = A0 - cB * pB
        return (4.0 * math.pi * D_B / psi) * (pB * R - gb) + V * (k_s * pA / psi - k_b) + Q * pA / psi

    R_scale = (3.0 * V_s * max(phibar, 1e-6) / (4.0 * math.pi * psi)) ** (1.0 / 3.0)
    R_lo_f = max(1e-4 * gb, 1e-6 * R_scale)
    R_hi_f = min(1.2 * R_scale, (3.0 * 0.99 * V_s / (4.0 * math.pi)) ** (1.0 / 3.0))
    gf = np.geomspace(R_lo_f, R_hi_f, 400)
    vf = [dV_full(float(r)) for r in gf]
    rf = [brentq(dV_full, float(gf[i]), float(gf[i + 1]), xtol=1e-9 * R_lo_f, rtol=1e-9) for i in range(len(gf) - 1) if vf[i] is not None and vf[i + 1] is not None and vf[i] * vf[i + 1] < 0.0]
    R_bar_full = (rf[-1] / w) if rf else None

    # Reaction-only Eq. (10) (J = 0, V << V_s)
    def dV_noJ(R):
        V = (4.0 * math.pi / 3.0) * (R ** 3)
        A0 = phibar - (V / V_s) * psi
        pA, pB = (k_b / k, A0 - k_b / k) if k * A0 > k_b else (A0, 0.0)
        return (4.0 * math.pi * D_B / psi) * (pB * R - gb) + V * (k_s * pA / psi - k_b) + Q * pA / psi

    R_max = max(500.0 * w, R_scale * 1.5)
    gn = np.geomspace(0.5 * w, R_max, 400)
    vn = [dV_noJ(float(r)) for r in gn]
    rn = [brentq(dV_noJ, float(gn[i]), float(gn[i + 1]), xtol=1e-10 * w, rtol=1e-10) / w for i in range(len(gn) - 1) if vn[i] * vn[i + 1] < 0.0]
    R_bar_noJ = rn[-1] if rn else None
    R_crit_fin = rn[0] if len(rn) >= 2 else None

    return {
        "R_crit_inf_norm": R_crit_inf,
        "R_bar_eq23_norm": R_bar_eq23,
        "R_bar_quad_appB_norm": R_bar_quad,
        "R_bar_eqB1_norm": R_bar_eqB1,
        "R_bar_full_norm": R_bar_full,
        "R_bar_noJ_norm": R_bar_noJ,
        "R_crit_fin_norm": R_crit_fin,
    }
