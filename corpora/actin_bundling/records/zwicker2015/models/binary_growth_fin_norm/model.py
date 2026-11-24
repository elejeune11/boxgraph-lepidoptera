import math


def run(inputs, settings):
    R_norm = inputs["R_norm"]
    kf_n = inputs["k_f_norm"]
    kb_n = inputs["k_b_norm"]
    phi_B_bar = inputs["phi_B_bar"]
    Vs_norm = inputs["Vs_norm"]
    w = inputs["w"]
    D_B = inputs["D_B"]

    k_0 = D_B / (w ** 2)
    k_f = kf_n * k_0
    k_b = kb_n * k_0
    gb = w / 6.0
    R = R_norm * w
    V_s = Vs_norm * (w ** 3)
    V = (4.0 * math.pi / 3.0) * (R ** 3)
    if kf_n == 0.0 and kb_n == 0.0:
        if V >= V_s:
            return {"dR_dt_norm": None}
        phi_0_B = (V_s * phi_B_bar - V) / (V_s - V)
    else:
        phi_0_B = (k_f + 4.0 * math.pi * D_B * gb / V_s) / (k_f + k_b + 4.0 * math.pi * D_B * R / V_s)
    dR_dt = D_B * (phi_0_B / R - gb / (R ** 2)) - k_b * R / 3.0
    return {"dR_dt_norm": dR_dt / (w * k_0)}
