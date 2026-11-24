import math


def run(inputs, settings):
    k_f = inputs["k_f"]
    k_b = inputs["k_b"]
    D_B = inputs["D_B"]
    gb = inputs["gamma_beta"]
    R_bar = inputs["R_bar"]
    phi_0_B = inputs["phi_0_B"]
    V_s = inputs["V_s"]
    N = inputs["N"]

    J00 = -k_f - k_b - 4.0 * math.pi * D_B * N * R_bar / V_s
    # Printed Eq. (D3b) includes a factor of N here
    J0i = -4.0 * math.pi * D_B * N * phi_0_B / V_s
    Ji0 = D_B / R_bar
    Jii = D_B * (2.0 * gb / R_bar**3 - phi_0_B / R_bar**2) - k_b / 3.0

    disc = (J00 - Jii) ** 2 + 4.0 * N * J0i * Ji0
    if disc < 0.0:
        return {"e_plus": None, "e_minus": None, "e_tilde": Jii}
    sq = math.sqrt(disc)
    return {
        "e_plus": 0.5 * (J00 + Jii + sq),
        "e_minus": 0.5 * (J00 + Jii - sq),
        "e_tilde": Jii,
    }


def symbolic():
    import sympy as sp

    k_f, k_b, D_B, gamma_beta, R_bar, phi_0_B, V_s, N = sp.symbols(
        "k_f k_b D_B gamma_beta R_bar phi_0_B V_s N", positive=True
    )
    J00 = -k_f - k_b - 4 * sp.pi * D_B * N * R_bar / V_s
    J0i = -4 * sp.pi * D_B * N * phi_0_B / V_s
    Ji0 = D_B / R_bar
    Jii = D_B * (2 * gamma_beta / R_bar**3 - phi_0_B / R_bar**2) - k_b / 3
    sq = sp.sqrt((J00 - Jii) ** 2 + 4 * N * J0i * Ji0)
    return {
        "e_plus": sp.Rational(1, 2) * (J00 + Jii + sq),
        "e_minus": sp.Rational(1, 2) * (J00 + Jii - sq),
        "e_tilde": Jii,
    }
