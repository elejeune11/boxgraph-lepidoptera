def run(inputs, settings):
    w = inputs["w"]
    D_B = inputs["D_B"]
    psi = inputs["psi"]
    k_0 = D_B / (w**2)
    k_f = inputs["kf_norm"] * k_0
    k_b = inputs["kb_norm"] * k_0
    k = inputs["k_norm"] * k_0
    return {
        "k_0": k_0,
        "k_f": k_f,
        "k_b": k_b,
        "k": k,
        "k_s": k_f + k * psi,
        "V_s": inputs["Vs_norm"] * (w**3),
        "Q": inputs["Q_norm"] * (w**3) * k_0,
        "gamma_beta_continuum": w * psi / 6.0,
        "gamma_beta_capillary": w / 6.0,
    }


def symbolic():
    import sympy as sp

    kf_norm, kb_norm, k_norm, Vs_norm, Q_norm, psi, w, D_B = sp.symbols(
        "kf_norm kb_norm k_norm Vs_norm Q_norm psi w D_B", positive=True
    )
    k_0 = D_B / w**2
    k_f = kf_norm * k_0
    k = k_norm * k_0
    return {
        "k_0": k_0,
        "k_f": k_f,
        "k_b": kb_norm * k_0,
        "k": k,
        "k_s": k_f + k * psi,
        "V_s": Vs_norm * w**3,
        "Q": Q_norm * w**3 * k_0,
        "gamma_beta_continuum": w * psi / 6,
        "gamma_beta_capillary": w / 6,
    }
