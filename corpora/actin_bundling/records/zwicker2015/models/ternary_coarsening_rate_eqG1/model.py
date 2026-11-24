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

    R_bar = (3.0 * V_bar / (4.0 * math.pi)) ** (1.0 / 3.0)
    # Eq. (G3)
    phi_0_B = gb / R_bar - J / (4.0 * math.pi * R_bar * D_B)
    # Eq. (G2)
    phi_0_A = (
        (k_b * psi * V_bar + 4.0 * math.pi * D_B * (gb - phi_0_B * R_bar))
        / (k_s * V_bar)
        * (1.0 - Q / (k_s * V_bar))
    )
    # Eq. (G1)
    lam = D_B * phi_0_B / (psi * R_bar**2) + k_s * phi_0_A / psi - k_b
    return {"lam": lam, "phi_0_A": phi_0_A, "phi_0_B": phi_0_B}


def symbolic():
    import sympy as sp

    D_B, gamma_beta, psi, V_bar, J, Q, k_b, k_s = sp.symbols(
        "D_B gamma_beta psi V_bar J Q k_b k_s", positive=True
    )
    R_bar = (3 * V_bar / (4 * sp.pi)) ** sp.Rational(1, 3)
    phi_0_B = gamma_beta / R_bar - J / (4 * sp.pi * R_bar * D_B)
    phi_0_A = (
        (k_b * psi * V_bar + 4 * sp.pi * D_B * (gamma_beta - phi_0_B * R_bar))
        / (k_s * V_bar)
        * (1 - Q / (k_s * V_bar))
    )
    lam = sp.simplify(D_B * phi_0_B / (psi * R_bar**2) + k_s * phi_0_A / psi - k_b)
    return {"lam": lam, "phi_0_A": sp.simplify(phi_0_A), "phi_0_B": phi_0_B}
