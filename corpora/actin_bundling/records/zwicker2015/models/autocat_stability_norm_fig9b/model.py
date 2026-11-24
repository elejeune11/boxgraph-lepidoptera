import math
import numpy as np
from scipy.optimize import brentq


def run(inputs, settings):
    phibar = inputs["phibar"]
    k_n = inputs["k_norm"]
    kb_n = inputs["k_b_norm"]
    psi = inputs["psi"]
    Vs_n = inputs["Vs_norm"]
    w = inputs["w"]
    D_B = inputs["D_B"]

    k_0 = D_B / (w ** 2)
    k = k_n * k_0
    k_b = kb_n * k_0
    k_s = k * psi
    V_s = Vs_n * (w ** 3)
    gb_cap = w / 6.0
    gb_cont = w * psi / 6.0
    Q_unit = (w ** 3) * k_0

    Qstab_cap = (4.0 * math.pi * D_B * k_s * gb_cap / (psi * k_b)) / Q_unit
    Qstab_cont = (4.0 * math.pi * D_B * k_s * gb_cont / (psi * k_b)) / Q_unit
    Qstab_eq24_printed = (4.0 * math.pi * D_B * k * gb_cap / (psi * k_b)) / Q_unit

    def bg_full(r1, r2):
        Vtot = (4.0 * math.pi / 3.0) * (r1 ** 3 + r2 ** 3)
        A0 = phibar - (Vtot / V_s) * psi
        cB = (V_s - Vtot) / V_s
        n_act = (1 if r1 > 0.0 else 0) + (1 if r2 > 0.0 else 0)
        Rsum = (r1 if r1 > 0.0 else 0.0) + (r2 if r2 > 0.0 else 0.0)
        qa = -k * cB
        qb = -k_b + k * A0 - 4.0 * math.pi * D_B * Rsum / (V_s - Vtot)
        qc = 4.0 * math.pi * D_B * n_act * gb_cap / (V_s - Vtot)
        disc = qb ** 2 - 4.0 * qa * qc
        pB = max(0.0, (-qb - math.sqrt(max(disc, 0.0))) / (2.0 * qa))
        return A0 - cB * pB, pB

    def bg_no_kAB(r1, r2):
        Vtot = (4.0 * math.pi / 3.0) * (r1 ** 3 + r2 ** 3)
        n_act = (1 if r1 > 0.0 else 0) + (1 if r2 > 0.0 else 0)
        Rsum = (r1 if r1 > 0.0 else 0.0) + (r2 if r2 > 0.0 else 0.0)
        pB = (4.0 * math.pi * D_B * n_act * gb_cap / (V_s - Vtot)) / (k_b + 4.0 * math.pi * D_B * Rsum / (V_s - Vtot))
        pA = phibar - ((V_s - Vtot) / V_s) * pB - (Vtot / V_s) * psi
        return pA, pB

    def bg_noJ(r1, r2):
        Vtot = (4.0 * math.pi / 3.0) * (r1 ** 3 + r2 ** 3)
        A0 = phibar - (Vtot / V_s) * psi
        return (k_b / k, A0 - k_b / k) if k * A0 > k_b else (A0, 0.0)

    def find_Q_multi(bg_fn):
        if phibar <= k_b / k:
            return None

        def lam_of_Q(q_n):
            Q = q_n * Q_unit
            Rs = np.geomspace(1.0 * w, 200.0 * w, 200)

            def dV_sym(r):
                pA, pB = bg_fn(r, r)
                V = (4.0 * math.pi / 3.0) * (r ** 3)
                return (4.0 * math.pi * D_B / psi) * (pB * r - gb_cap) + V * (k_s * pA / psi - k_b) + Q * pA / psi

            vals = [dV_sym(float(r)) for r in Rs]
            Rbar = None
            for i in range(len(Rs) - 1, 0, -1):
                if vals[i - 1] * vals[i] <= 0.0:
                    Rbar = brentq(dV_sym, float(Rs[i - 1]), float(Rs[i]), xtol=1e-9 * w, rtol=1e-9)
                    break
            if Rbar is None:
                return 1.0
            h = Rbar * 1e-5
            pAp, pBp = bg_fn(Rbar + h, Rbar - h)
            Vp = (4.0 * math.pi / 3.0) * ((Rbar + h) ** 3)
            Vm = (4.0 * math.pi / 3.0) * ((Rbar - h) ** 3)
            dv1 = (4.0 * math.pi * D_B / psi) * (pBp * (Rbar + h) - gb_cap) + Vp * (k_s * pAp / psi - k_b) + Q * pAp / psi
            dv2 = (4.0 * math.pi * D_B / psi) * (pBp * (Rbar - h) - gb_cap) + Vm * (k_s * pAp / psi - k_b) + Q * pAp / psi
            return (dv1 / (4.0 * math.pi * (Rbar + h) ** 2) - dv2 / (4.0 * math.pi * (Rbar - h) ** 2)) / (2.0 * h)

        qgrid = np.linspace(40.0, 260.0, 23)
        lvals = [lam_of_Q(float(q)) for q in qgrid]
        for i in range(len(qgrid) - 1):
            if lvals[i] * lvals[i + 1] <= 0.0:
                return float(brentq(lam_of_Q, float(qgrid[i]), float(qgrid[i + 1]), xtol=1e-4, rtol=1e-6))
        return None

    def find_Q_droplet(bg_fn):
        Rs = np.geomspace(1.0 * w, 180.0 * w, 300)
        qmin = None
        for r in Rs:
            r = float(r)
            pA, pB = bg_fn(r, 0.0)
            if pA <= 1e-12:
                continue
            V = (4.0 * math.pi / 3.0) * (r ** 3)
            q = (-(4.0 * math.pi * D_B / psi) * (pB * r - gb_cap) - V * (k_s * pA / psi - k_b)) * psi / (pA * Q_unit)
            if qmin is None or q < qmin:
                qmin = q
        return max(0.0, float(qmin)) if qmin is not None else None

    return {
        "Qstab_cap_norm": Qstab_cap,
        "Qstab_cont_norm": Qstab_cont,
        "Qstab_eq24_printed_norm": Qstab_eq24_printed,
        "Q_multi_full_norm": find_Q_multi(bg_full),
        "Q_multi_no_kAB_norm": find_Q_multi(bg_no_kAB),
        "Q_droplet_full_norm": find_Q_droplet(bg_full),
        "Q_droplet_noJ_norm": find_Q_droplet(bg_noJ),
    }
