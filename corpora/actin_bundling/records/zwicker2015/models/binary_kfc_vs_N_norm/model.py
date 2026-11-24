def run(inputs, settings):
    _ = inputs["N"]
    ratio = inputs["kb_over_kf"]
    w = inputs["w"]
    D_B = inputs["D_B"]
    k_0 = D_B / (w ** 2)
    gb = w / 6.0
    phi_0_B = 1.0 / (1.0 + ratio)
    kfc = (4.0 / 9.0) * (D_B / (gb ** 2)) * (phi_0_B ** 4) / (1.0 - phi_0_B)
    return {"k_f_c_norm": kfc / k_0}
