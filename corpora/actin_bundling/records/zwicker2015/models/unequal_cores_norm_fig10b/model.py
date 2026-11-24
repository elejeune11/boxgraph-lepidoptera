import math
import numpy as np
from scipy.optimize import brentq, root


def run(inputs, settings):
    Q1_n = inputs["Q1_norm"]
    Q2_n = inputs["Q2_norm"]
    k_n = inputs["k_norm"]
    kb_n = inputs["k_b_norm"]
    psi = inputs["psi"]
    phibar = inputs["phibar"]
    Vs_n = inputs["Vs_norm"]
    w = inputs["w"]
    D_B = inputs["D_B"]

    k_0 = D_B / (w ** 2)
    k = k_n * k_0
    k_b = kb_n * k_0
    k_s = k * psi
    gb = w / 6.0
    V_s = Vs_n * (w ** 3)
    Q_unit = (w ** 3) * k_0
    Q1 = Q1_n * Q_unit
    Q2 = Q2_n * Q_unit

    def dV_full(r1, r2):
        Vtot = (4.0 * math.pi / 3.0) * (r1 ** 3 + r2 ** 3)
        if Vtot >= 0.999999 * V_s:
            return None, None
        A0 = phibar - (Vtot / V_s) * psi
        cB = (V_s - Vtot) / V_s
        n_act = (1 if r1 > 0.0 else 0) + (1 if r2 > 0.0 else 0)
        Rsum = (r1 if r1 > 0.0 else 0.0) + (r2 if r2 > 0.0 else 0.0)
        qa = -k * cB
        qb = -k_b + k * A0 - 4.0 * math.pi * D_B * Rsum / (V_s - Vtot)
        qc = 4.0 * math.pi * D_B * n_act * gb / (V_s - Vtot)
        disc = qb ** 2 - 4.0 * qa * qc
        if disc < 0.0:
            return None, None
        pB = max(0.0, (-qb - math.sqrt(disc)) / (2.0 * qa))
        pA = A0 - cB * pB
        v1 = (4.0 * math.pi / 3.0) * (r1 ** 3)
        v2 = (4.0 * math.pi / 3.0) * (r2 ** 3)
        dv1 = (4.0 * math.pi * D_B / psi) * (pB * r1 - gb) + v1 * (k_s * pA / psi - k_b) + Q1 * pA / psi
        dv2 = (4.0 * math.pi * D_B / psi) * (pB * r2 - gb) + v2 * (k_s * pA / psi - k_b) + Q2 * pA / psi
        return dv1, dv2

    # 1. Full Eq. (10)+(12)
    V_est = max((V_s / psi) * max(phibar - k_b / k, 1e-9), 1e-30)
    R_sc = (3.0 * V_est / (4.0 * math.pi)) ** (1.0 / 3.0)
    rate_sc = max(k_b * R_sc, 1e-30)
    stable_pair = False
    R1s = R2s = None
    if abs(Q1 - Q2) <= 1e-12 * max(abs(Q1), abs(Q2), 1.0):
        g = np.geomspace(1.0 * w, 200.0 * w, 300)
        v = [dV_full(float(r), float(r))[0] for r in g]
        rts = [brentq(lambda r: dV_full(r, r)[0], float(g[i]), float(g[i + 1]), xtol=1e-9 * w, rtol=1e-9) for i in range(len(g) - 1) if v[i] is not None and v[i + 1] is not None and v[i] * v[i + 1] < 0.0]
        if rts:
            R1s = R2s = rts[-1]
            h = R1s * 1e-5
            dp1, dp2 = dV_full(R1s + h, R1s - h)
            lam_a = (dp1 / (4.0 * math.pi * (R1s + h) ** 2) - dp2 / (4.0 * math.pi * (R1s - h) ** 2)) / (2.0 * h)
            ds1, _ = dV_full(R1s + h, R1s + h)
            stable_pair = (lam_a < 0.0 and ds1 < 0.0)
    else:
        f1 = Q1 / (Q1 + Q2)
        r1g = (3.0 * V_est * f1 / (4.0 * math.pi)) ** (1.0 / 3.0)
        r2g = (3.0 * V_est * (1.0 - f1) / (4.0 * math.pi)) ** (1.0 / 3.0)

        def res_u(u):
            r1 = R_sc * math.exp(max(min(float(u[0]), 20.0), -20.0))
            r2 = R_sc * math.exp(max(min(float(u[1]), 20.0), -20.0))
            dv1, dv2 = dV_full(r1, r2)
            if dv1 is None or dv2 is None:
                return [-1e6, -1e6]
            return [dv1 / (4.0 * math.pi * r1 ** 2 * rate_sc), dv2 / (4.0 * math.pi * r2 ** 2 * rate_sc)]

        sol = root(res_u, [math.log(r1g / R_sc), math.log(r2g / R_sc)], tol=1e-9)
        if sol.success:
            R1s = R_sc * math.exp(float(sol.x[0]))
            R2s = R_sc * math.exp(float(sol.x[1]))
            h = 1e-5 * min(R1s, R2s)
            J2 = np.zeros((2, 2))
            for j, (dr1, dr2) in enumerate([(h, 0.0), (0.0, h)]):
                dp1, dp2 = dV_full(R1s + dr1, R2s + dr2)
                dm1, dm2 = dV_full(R1s - dr1, R2s - dr2)
                J2[0, j] = (dp1 / (4.0 * math.pi * (R1s + dr1) ** 2) - dm1 / (4.0 * math.pi * (R1s - dr1) ** 2)) / (2.0 * h)
                J2[1, j] = (dp2 / (4.0 * math.pi * (R2s + dr2) ** 2) - dm2 / (4.0 * math.pi * (R2s - dr2) ** 2)) / (2.0 * h)
            stable_pair = float(np.max(np.linalg.eigvals(J2).real)) < 0.0

    if stable_pair:
        V1_full = (4.0 * math.pi / 3.0) * ((R1s / w) ** 3)
        V2_full = (4.0 * math.pi / 3.0) * ((R2s / w) ** 3)
    else:
        g = np.geomspace(1.0 * w, 200.0 * w, 300)
        v = [dV_full(0.0, float(r))[1] for r in g]
        rts = [brentq(lambda r: dV_full(0.0, r)[1], float(g[i]), float(g[i + 1]), xtol=1e-9 * w, rtol=1e-9) for i in range(len(g) - 1) if v[i] is not None and v[i + 1] is not None and v[i] * v[i + 1] < 0.0]
        V1_full = 0.0
        V2_full = (4.0 * math.pi / 3.0) * ((rts[-1] / w) ** 3) if rts else 0.0
    Vtot_full = V1_full + V2_full

    # 2. Reaction-only Eq. (10) (exact closed-form quadratic!)
    Q_stab = 4.0 * math.pi * D_B * k * gb / k_b
    if Q1 <= Q_stab:
        A1 = k * V_s
        B1 = -(V_s * (k_b + k * phibar) + Q2)
        C1 = V_s * phibar * k_b + 4.0 * math.pi * D_B * gb
        pA1 = (-B1 - math.sqrt(max(B1 ** 2 - 4.0 * A1 * C1, 0.0))) / (2.0 * A1)
        V1_noJ = 0.0
        V2_noJ = (Q2 * pA1 - 4.0 * math.pi * D_B * gb) / (psi * (k_b - k * pA1)) / (w ** 3)
    else:
        A2 = k * V_s
        B2 = -(V_s * (k_b + k * phibar) + Q1 + Q2)
        C2 = V_s * phibar * k_b + 8.0 * math.pi * D_B * gb
        pA2 = (-B2 - math.sqrt(max(B2 ** 2 - 4.0 * A2 * C2, 0.0))) / (2.0 * A2)
        V1_noJ = (Q1 * pA2 - 4.0 * math.pi * D_B * gb) / (psi * (k_b - k * pA2)) / (w ** 3)
        V2_noJ = (Q2 * pA2 - 4.0 * math.pi * D_B * gb) / (psi * (k_b - k * pA2)) / (w ** 3)
    Vtot_noJ = V1_noJ + V2_noJ

    # 3. Zero-surface-tension limit (V1/V2 = Q1/Q2 with Eq. 23)
    Vtot_eq23 = ((V_s / psi) * (phibar - k_b / k)) / (w ** 3)
    ratio = Q1 / Q2
    V1_zg = Vtot_eq23 * ratio / (1.0 + ratio)
    V2_zg = Vtot_eq23 / (1.0 + ratio)

    return {
        "V1_full_norm": V1_full,
        "V2_full_norm": V2_full,
        "Vtot_full_norm": Vtot_full,
        "V1_noJ_norm": V1_noJ,
        "V2_noJ_norm": V2_noJ,
        "Vtot_noJ_norm": Vtot_noJ,
        "V1_zerogamma_norm": V1_zg,
        "V2_zerogamma_norm": V2_zg,
        "Vtot_eq23_norm": Vtot_eq23,
    }
