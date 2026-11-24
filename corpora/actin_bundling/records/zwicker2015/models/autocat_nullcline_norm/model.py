import math
import numpy as np
from scipy.optimize import brentq, root


def run(inputs, settings):
    R1_n = inputs["R1_norm"]
    k_n = inputs["k_norm"]
    kb_n = inputs["k_b_norm"]
    psi = inputs["psi"]
    phibar = inputs["phibar"]
    Vs_n = inputs["Vs_norm"]
    Q1_n = inputs["Q1_norm"]
    Q2_n = inputs["Q2_norm"]
    w = inputs["w"]
    D_B = inputs["D_B"]

    k_0 = D_B / (w ** 2)
    k = k_n * k_0
    k_b = kb_n * k_0
    k_s = k * psi
    gb = w / 6.0
    V_s = Vs_n * (w ** 3)
    Q1 = Q1_n * (w ** 3) * k_0
    Q2 = Q2_n * (w ** 3) * k_0
    R1 = R1_n * w

    def bg_full(r1, r2):
        Vtot = (4.0 * math.pi / 3.0) * (r1 ** 3 + r2 ** 3)
        if Vtot >= 0.9999 * V_s:
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
        return A0 - cB * pB, pB

    def bg_noJ(r1, r2):
        Vtot = (4.0 * math.pi / 3.0) * (r1 ** 3 + r2 ** 3)
        A0 = phibar - (Vtot / V_s) * psi
        return (k_b / k, A0 - k_b / k) if k * A0 > k_b else (A0, 0.0)

    def dV_pair(r1, r2, bg):
        pA, pB = bg(r1, r2)
        if pA is None:
            return None, None
        v1 = (4.0 * math.pi / 3.0) * (r1 ** 3)
        v2 = (4.0 * math.pi / 3.0) * (r2 ** 3)
        dv1 = (4.0 * math.pi * D_B / psi) * (pB * r1 - gb) + v1 * (k_s * pA / psi - k_b) + Q1 * pA / psi
        dv2 = (4.0 * math.pi * D_B / psi) * (pB * r2 - gb) + v2 * (k_s * pA / psi - k_b) + Q2 * pA / psi
        return dv1, dv2

    def find_R2(bg, idx):
        grid = np.linspace(0.5 * w, 200.0 * w, 250)
        vals = [dV_pair(R1, float(r2), bg)[idx] for r2 in grid]
        rts = []
        for i in range(len(grid) - 1):
            a, b = vals[i], vals[i + 1]
            if a is not None and b is not None and a * b <= 0.0:
                if a == 0.0:
                    rts.append(float(grid[i]) / w)
                elif b != 0.0:
                    rts.append(float(brentq(lambda r2: dV_pair(R1, r2, bg)[idx], float(grid[i]), float(grid[i + 1]), xtol=1e-9 * w, rtol=1e-9)) / w)
        return rts[-1] if rts else None

    # Fixed points (symmetric vs unequal)
    q_scale = max(abs(Q1), abs(Q2))
    is_sym = (q_scale == 0.0) or (abs(Q1 - Q2) <= 1e-12 * q_scale)

    R2_fp_full = None
    R2_fp_noJ = None
    if is_sym:
        g2 = np.geomspace(1.0 * w, 200.0 * w, 300)
        vf = [dV_pair(float(r), float(r), bg_full)[0] for r in g2]
        rf = [brentq(lambda r: dV_pair(r, r, bg_full)[0], float(g2[i]), float(g2[i + 1]), xtol=1e-10 * w, rtol=1e-10) / w for i in range(len(g2) - 1) if vf[i] is not None and vf[i + 1] is not None and vf[i] * vf[i + 1] < 0.0]
        R2_fp_full = rf[-1] if rf else None
        vn = [dV_pair(float(r), float(r), bg_noJ)[0] for r in g2]
        rn = [brentq(lambda r: dV_pair(r, r, bg_noJ)[0], float(g2[i]), float(g2[i + 1]), xtol=1e-10 * w, rtol=1e-10) / w for i in range(len(g2) - 1) if vn[i] * vn[i + 1] < 0.0]
        R2_fp_noJ = rn[-1] if rn else None
    else:
        V_est = max((V_s / psi) * max(phibar - k_b / k, 1e-9), 1e-30)
        f1 = Q1 / (Q1 + Q2)
        r1g = (3.0 * V_est * f1 / (4.0 * math.pi)) ** (1.0 / 3.0)
        r2g = (3.0 * V_est * (1.0 - f1) / (4.0 * math.pi)) ** (1.0 / 3.0)
        R_sc = (3.0 * V_est / (4.0 * math.pi)) ** (1.0 / 3.0)
        rate_sc = max(k_b * R_sc, 1e-30)

        def res_f(u):
            r1 = R_sc * math.exp(max(min(float(u[0]), 15.0), -15.0))
            r2 = R_sc * math.exp(max(min(float(u[1]), 15.0), -15.0))
            dv1, dv2 = dV_pair(r1, r2, bg_full)
            if dv1 is None or dv2 is None:
                return [-1e6, -1e6]
            return [dv1 / (4.0 * math.pi * r1 ** 2 * rate_sc), dv2 / (4.0 * math.pi * r2 ** 2 * rate_sc)]

        sol = root(res_f, [math.log(r1g / R_sc), math.log(r2g / R_sc)], tol=1e-9)
        if sol.success:
            R2_fp_full = (R_sc * math.exp(float(sol.x[1]))) / w

        # Exact closed-form two-droplet steady state for Q1, Q2 > Q_stab (phi_0^B = 0)
        A = k * V_s
        B = -(V_s * (k_b + k * phibar) + Q1 + Q2)
        C = V_s * phibar * k_b + 8.0 * math.pi * D_B * gb
        disc = B ** 2 - 4.0 * A * C
        if disc >= 0.0:
            pA = (-B - math.sqrt(disc)) / (2.0 * A)
            if pA < k_b / k and Q1 * pA > 4.0 * math.pi * D_B * gb and Q2 * pA > 4.0 * math.pi * D_B * gb:
                V2 = (Q2 * pA - 4.0 * math.pi * D_B * gb) / (psi * (k_b - k * pA))
                R2_fp_noJ = ((3.0 * V2 / (4.0 * math.pi)) ** (1.0 / 3.0)) / w

    return {
        "R2_null1_full_norm": find_R2(bg_full, 0),
        "R2_null2_full_norm": find_R2(bg_full, 1),
        "R2_fp_full_norm": R2_fp_full,
        "R2_null1_noJ_norm": find_R2(bg_noJ, 0),
        "R2_null2_noJ_norm": find_R2(bg_noJ, 1),
        "R2_fp_noJ_norm": R2_fp_noJ,
    }
