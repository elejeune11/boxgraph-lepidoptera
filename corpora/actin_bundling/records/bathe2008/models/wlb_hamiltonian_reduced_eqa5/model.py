"""Appendix Eq. (A5): effective bending stiffness kappa_B = 2 * H_j[r_j, Delta_u*] / (q_j^4 * r_j^2)."""


def run(inputs, settings):
    M = inputs["M"]
    kappa_f = inputs["kappa_f"]
    E_f_A_f = inputs["E_f_A_f"]
    d_f = inputs["d_f"]
    t = inputs["t"]
    k_x = inputs["k_x"]
    delta = inputs["delta"]
    L = inputs["L"]
    q_j = inputs["q_j"]
    N = 4.0 * M * M
    alpha = k_x * (L ** 2) / (E_f_A_f * delta)
    # From Eq. (A5): H_j / (q_j^4 * r_j^2) = (1/2)*N*kappa_f + M^2*E_f*A_f*(d_f+t)^2*(4M^2-1) / (6 + (q_j*L)^2*M*(2M+1)/alpha)
    # Since H_j = (1/2) * kappa_B * q_j^4 * r_j^2, multiply by 2:
    kappa_B = N * kappa_f + 2.0 * (M ** 2) * E_f_A_f * ((d_f + t) ** 2) * (4.0 * M * M - 1.0) / (
        6.0 + ((q_j * L) ** 2) * M * (2.0 * M + 1.0) / alpha
    )
    return {"kappa_B": kappa_B}


def symbolic():
    import sympy as sp

    M, kappa_f, E_f_A_f, d_f, t, k_x, delta, L, q_j = sp.symbols(
        "M kappa_f E_f_A_f d_f t k_x delta L q_j", positive=True
    )
    N = 4 * M**2
    alpha = k_x * L**2 / (E_f_A_f * delta)
    kappa_B = N * kappa_f + 2 * M**2 * E_f_A_f * (d_f + t) ** 2 * (4 * M**2 - 1) / (
        6 + (q_j * L) ** 2 * M * (2 * M + 1) / alpha
    )
    return {"kappa_B": kappa_B}
