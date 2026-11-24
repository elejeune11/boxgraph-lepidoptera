import math


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
        dV_c = (Q * k_b - 4.0 * math.pi * D_B * k * gb) / ((k * phibar - k_b) * k * psi)
        if V0 + dV_c > 0.0:
            R_bar_eqB1 = ((3.0 * (V0 + dV_c) / (4.0 * math.pi)) ** (1.0 / 3.0)) / w

    A = k * (psi ** 2)
    B = -V_s * psi * (k * phibar - k_b - Q / V_s)
    C = -V_s * (Q * phibar - 4.0 * math.pi * D_B * gb)
    disc = B ** 2 - 4.0 * A * C
    R_bar_quad = ((3.0 * ((-B + math.sqrt(disc)) / (2.0 * A)) / (4.0 * math.pi)) ** (1.0 / 3.0)) / w if (disc >= 0.0 and (-B + math.sqrt(disc)) > 0.0) else None

    return {
        "R_crit_inf_norm": R_crit_inf,
        "R_bar_eq23_norm": R_bar_eq23,
        "R_bar_quad_appB_norm": R_bar_quad,
        "R_bar_eqB1_norm": R_bar_eqB1,
        "R_bar_full_norm": None,
        "R_bar_noJ_norm": R_bar_quad,
        "R_crit_fin_norm": None,
    }
