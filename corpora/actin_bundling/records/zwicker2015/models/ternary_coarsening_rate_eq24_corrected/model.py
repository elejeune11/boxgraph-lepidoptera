import math


def run(inputs, settings):
    D_B = inputs["D_B"]
    gb = inputs["gamma_beta"]
    psi = inputs["psi"]
    V_bar = inputs["V_bar"]
    J = inputs["J"]
    Q = inputs["Q"]
    k_b = inputs["k_b"]
    k_s = inputs["k_s"]

    lam = (
        4.0 * math.pi * D_B * gb / (3.0 * psi * V_bar)
        + 2.0 * J / (3.0 * psi * V_bar)
        - (Q / (V_bar * k_s)) * (k_b + J / (psi * V_bar))
    )
    return {"lam": lam}


def symbolic():
    import sympy as sp

    D_B, gamma_beta, psi, V_bar, J, Q, k_b, k_s = sp.symbols(
        "D_B gamma_beta psi V_bar J Q k_b k_s", positive=True
    )
    lam = (
        4 * sp.pi * D_B * gamma_beta / (3 * psi * V_bar)
        + 2 * J / (3 * psi * V_bar)
        - (Q / (V_bar * k_s)) * (k_b + J / (psi * V_bar))
    )
    return {"lam": lam}
