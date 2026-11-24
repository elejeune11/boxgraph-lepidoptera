import math


def run(inputs, settings):
    D_B = inputs["D_B"]
    k_s = inputs["k_s"]
    R = inputs["R"]
    gb = inputs["gamma_beta"]
    phi_0_B = inputs["phi_0_B"]
    psi = inputs["psi"]
    k_b = inputs["k_b"]

    phi_plus_B = gb / R
    Q_stab = (4.0 * math.pi * D_B * k_s * R * (3.0 * phi_plus_B - 2.0 * phi_0_B)) / (3.0 * psi * k_b)
    Q_stab_inf = (4.0 * math.pi * D_B * k_s * gb) / (psi * k_b)
    return {"Q_stab": Q_stab, "Q_stab_inf": Q_stab_inf}


def symbolic():
    import sympy as sp

    D_B, k_s, R, gamma_beta, phi_0_B, psi, k_b = sp.symbols(
        "D_B k_s R gamma_beta phi_0_B psi k_b", positive=True
    )
    phi_plus_B = gamma_beta / R
    Q_stab = 4 * sp.pi * D_B * k_s * R * (3 * phi_plus_B - 2 * phi_0_B) / (3 * psi * k_b)
    Q_stab_inf = 4 * sp.pi * D_B * k_s * gamma_beta / (psi * k_b)
    return {"Q_stab": sp.simplify(Q_stab), "Q_stab_inf": Q_stab_inf}
