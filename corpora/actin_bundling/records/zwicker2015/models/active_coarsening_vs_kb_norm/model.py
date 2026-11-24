import math


def run(inputs, settings):
    kb_n = inputs["k_b_norm"]
    ratio = inputs["kf_over_kb"]
    Vs_n = inputs["Vs_norm"]
    N = inputs["N"]
    w = inputs["w"]
    D_B = inputs["D_B"]
    k_0 = D_B / (w ** 2)
    k_b = kb_n * k_0
    k_f = ratio * k_b
    V_s = Vs_n * (w ** 3)
    V_bar = (k_f / (k_f + k_b)) * V_s / N
    R_bar = (3.0 * V_bar / (4.0 * math.pi)) ** (1.0 / 3.0)
    lam = D_B * (w / 6.0) / (R_bar ** 3) - (2.0 / 3.0) * k_b
    return {"lam_norm": lam / k_0}
