import math


def run(inputs, settings):
    V_bar_n = inputs["V_bar_norm"]
    w = inputs["w"]
    D_B = inputs["D_B"]
    k_0 = D_B / (w ** 2)
    V_bar = V_bar_n * (w ** 3)
    R_bar = (3.0 * V_bar / (4.0 * math.pi)) ** (1.0 / 3.0)
    lam = D_B * (w / 6.0) / (R_bar ** 3)
    return {"lam_norm": lam / k_0}
