def run(inputs, settings):
    D_B = inputs["D_B"]
    w = inputs["w"]
    k_f = inputs["k_f"]
    k_b = inputs["k_b"]
    q = inputs["q"]

    k_0 = D_B / (w**2)
    c_spin = (4.0 * k_f * k_b - k_f**2 - k_b**2) / ((k_f + k_b) ** 2)
    lam_q = 2.0 * D_B * q**2 * (c_spin - (w**2 / 8.0) * q**2) - (k_f + k_b)
    if c_spin > 0.0:
        lam_max = (4.0 * D_B / w**2) * (c_spin ** 2) - (k_f + k_b)
    else:
        lam_max = -(k_f + k_b)
    return {"lam_q": lam_q, "lam_q_norm": lam_q / k_0, "lam_max": lam_max}


def symbolic():
    import sympy as sp

    D_B, w, k_f, k_b, q = sp.symbols("D_B w k_f k_b q", positive=True)
    k_0 = D_B / w**2
    c_spin = (4 * k_f * k_b - k_f**2 - k_b**2) / (k_f + k_b) ** 2
    lam_q = 2 * D_B * q**2 * (c_spin - (w**2 / 8) * q**2) - k_f - k_b
    lam_max = (4 * D_B / w**2) * c_spin**2 - k_f - k_b
    return {"lam_q": lam_q, "lam_q_norm": sp.simplify(lam_q / k_0), "lam_max": lam_max}
