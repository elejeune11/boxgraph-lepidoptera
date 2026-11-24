def run(inputs, settings):
    w = inputs["w"]
    D_B = inputs["D_B"]
    k_0 = D_B / (w**2)
    return {
        "k_0": k_0,
        "k_f": inputs["kf_norm"] * k_0,
        "k_b": inputs["kb_norm"] * k_0,
        "V_s": inputs["Vs_norm"] * (w**3),
        "gamma_beta": w / 6.0,
    }


def symbolic():
    import sympy as sp

    kf_norm, kb_norm, Vs_norm, w, D_B = sp.symbols("kf_norm kb_norm Vs_norm w D_B", positive=True)
    k_0 = D_B / w**2
    return {
        "k_0": k_0,
        "k_f": kf_norm * k_0,
        "k_b": kb_norm * k_0,
        "V_s": Vs_norm * w**3,
        "gamma_beta": w / 6,
    }
