import math


def run(inputs, settings):
    R_n = inputs["R_norm"]
    k_n = inputs["k_norm"]
    kb_n = inputs["k_b_norm"]
    psi = inputs["psi"]
    phibar = inputs["phibar"]
    Vs_n = inputs["Vs_norm"]
    Q_n = inputs["Q_norm"]
    w = inputs["w"]
    D_B = inputs["D_B"]

    k_0 = D_B / (w ** 2)
    k = k_n * k_0
    k_b = kb_n * k_0
    k_s = k * psi
    gb = w / 6.0
    R = R_n * w
    V_s = Vs_n * (w ** 3)
    V = (4.0 * math.pi / 3.0) * (R ** 3)
    Q = Q_n * (w ** 3) * k_0

    def eval_dR(pA, pB):
        dV = (4.0 * math.pi * D_B / psi) * (pB * R - gb) + V * (k_s * pA / psi - k_b) + Q * pA / psi
        return dV / (4.0 * math.pi * (R ** 2) * w * k_0)

    pA_inf, pB_inf = (k_b / k, phibar - k_b / k) if k * phibar > k_b else (phibar, 0.0)
    dR_inf = eval_dR(pA_inf, pB_inf)

    dR_full = None
    if V < 0.999999 * V_s:
        A0 = phibar - (V / V_s) * psi
        cB = (V_s - V) / V_s
        qa = -k * cB
        qb = -k_b + k * A0 - 4.0 * math.pi * D_B * R / (V_s - V)
        qc = 4.0 * math.pi * D_B * gb / (V_s - V)
        disc = qb ** 2 - 4.0 * qa * qc
        if disc >= 0.0:
            pB_f = max(0.0, (-qb - math.sqrt(disc)) / (2.0 * qa))
            dR_full = eval_dR(A0 - cB * pB_f, pB_f)

    A0 = phibar - (V / V_s) * psi
    pA_n, pB_n = (k_b / k, A0 - k_b / k) if k * A0 > k_b else (A0, 0.0)
    dR_noJ = eval_dR(pA_n, pB_n)

    return {"dR_dt_inf_norm": dR_inf, "dR_dt_full_norm": dR_full, "dR_dt_noJ_norm": dR_noJ}
