def run(inputs, settings):
    R_norm = inputs["R_norm"]
    kf_n = inputs["k_f_norm"]
    kb_n = inputs["k_b_norm"]
    phi_B_bar = inputs["phi_B_bar"]
    w = inputs["w"]
    D_B = inputs["D_B"]

    k_0 = D_B / (w ** 2)
    k_f = kf_n * k_0
    k_b = kb_n * k_0
    gb = w / 6.0
    R = R_norm * w
    phi_0_B = phi_B_bar if (kf_n == 0.0 and kb_n == 0.0) else k_f / (k_f + k_b)
    dR_dt = D_B * (phi_0_B / R - gb / (R ** 2)) - k_b * R / 3.0
    return {"dR_dt_norm": dR_dt / (w * k_0)}
