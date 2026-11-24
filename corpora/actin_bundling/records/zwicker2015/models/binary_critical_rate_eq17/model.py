def run(inputs, settings):
    D_B = inputs["D_B"]
    gb = inputs["gamma_beta"]
    phi_0_B = inputs["phi_0_B"]
    k_f_c = (4.0 * D_B / (3.0 * gb) ** 2) * (phi_0_B ** 4) / (1.0 - phi_0_B)
    k_b_c = (4.0 * D_B / (3.0 * gb) ** 2) * (phi_0_B ** 3)
    return {"k_f_c": k_f_c, "k_b_c": k_b_c}


def symbolic():
    import sympy as sp

    D_B, gamma_beta, phi_0_B = sp.symbols("D_B gamma_beta phi_0_B", positive=True)
    k_f_c = (4 * D_B / (3 * gamma_beta) ** 2) * phi_0_B**4 / (1 - phi_0_B)
    k_b_c = (4 * D_B / (3 * gamma_beta) ** 2) * phi_0_B**3
    return {"k_f_c": k_f_c, "k_b_c": k_b_c}
