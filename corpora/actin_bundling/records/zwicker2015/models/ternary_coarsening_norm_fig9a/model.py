import math


def run(inputs, settings):
    Vb_n = inputs["V_bar_norm"]
    k_n = inputs["k_norm"]
    kb_n = inputs["k_b_norm"]
    psi = inputs["psi"]
    Q_n = inputs["Q_norm"]
    w = inputs["w"]
    D_B = inputs["D_B"]

    k_0 = D_B / (w ** 2)
    k_b = kb_n * k_0
    k_s = k_n * k_0 * psi
    Q = Q_n * (w ** 3) * k_0
    V_bar = Vb_n * (w ** 3)

    def eval_pair(gb):
        J_large = 4.0 * math.pi * D_B * gb
        lam_s = 4.0 * math.pi * D_B * gb / (3.0 * psi * V_bar) - (Q / (V_bar * k_s)) * k_b
        lam_l = 4.0 * math.pi * D_B * gb / (3.0 * psi * V_bar) + 2.0 * J_large / (3.0 * psi * V_bar) - (Q / (V_bar * k_s)) * (k_b + J_large / (psi * V_bar))
        return lam_s / k_0, lam_l / k_0

    ls_cont, ll_cont = eval_pair(w * psi / 6.0)
    ls_cap, ll_cap = eval_pair(w / 6.0)
    return {
        "lam_small_cont_norm": ls_cont,
        "lam_large_cont_norm": ll_cont,
        "lam_small_cap_norm": ls_cap,
        "lam_large_cap_norm": ll_cap,
    }
