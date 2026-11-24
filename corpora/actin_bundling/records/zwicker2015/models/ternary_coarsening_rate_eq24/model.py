import math


def run(inputs, settings):
    D_B = inputs["D_B"]
    gb = inputs["gamma_beta"]
    psi = inputs["psi"]
    V_bar = inputs["V_bar"]
    J = inputs["J"]
    Q = inputs["Q"]
    k_b = inputs["k_b"]
    k = inputs["k"]

    base = 4.0 * math.pi * D_B * gb / (3.0 * psi * V_bar) + 2.0 * J / (3.0 * psi * V_bar)
    if Q == 0.0:
        return {"lam": base}
    return {"lam": base - (Q / (V_bar * k)) * (k_b + J / (psi * V_bar))}


def symbolic():
    import sympy as sp

    D_B, gamma_beta, psi, V_bar, J, Q, k_b, k = sp.symbols(
        "D_B gamma_beta psi V_bar J Q k_b k", positive=True
    )
    lam = (
        4 * sp.pi * D_B * gamma_beta / (3 * psi * V_bar)
        + 2 * J / (3 * psi * V_bar)
        - (Q / (V_bar * k)) * (k_b + J / (psi * V_bar))
    )
    return {"lam": lam}
