import math


def run(inputs, settings):
    alpha = inputs["alpha"]
    xi_f0 = inputs["xi_f0_over_gamma"]
    D_tilde = inputs["D_tilde"]
    k_turn = inputs["k_tilde"]
    nu_bar = inputs["nu_bar"]
    eta = inputs["eta_s_over_eta_b"]

    # Eq. (S84): k_tilde_n^c
    num_c = alpha - k_turn + nu_bar * (1.0 + xi_f0)
    den_c = D_tilde + nu_bar * (1.0 + xi_f0)
    k_c = math.sqrt(num_c / den_c) if (num_c > 0.0 and den_c > 0.0) else 0.0

    # Eq. (S85): k_tilde_n^(1) = 1
    k_1 = 1.0

    # Eq. (S86): k_tilde_n^(2)
    num_2 = alpha - (1.0 + eta) * k_turn
    den_2 = (1.0 + eta) * D_tilde
    k_2 = math.sqrt(num_2 / den_2) if (num_2 > 0.0 and den_2 > 0.0) else 0.0

    # Unstable if Tr > 0 (k < k_c) or det < 0 (k between min(1, k_2) and max(1, k_2))
    k_crit = max(k_1, k_c, k_2)
    return {
        "k_tilde_trace_zero": k_c,
        "k_tilde_det_root1": k_1,
        "k_tilde_det_root2": k_2,
        "L_crit_over_2pi_r0": 1.0 / k_crit,
    }


def symbolic():
    import sympy as sp

    alpha, xi_f0_over_gamma, D_tilde, k_tilde, nu_bar, eta_s_over_eta_b = sp.symbols(
        "alpha xi_f0_over_gamma D_tilde k_tilde nu_bar eta_s_over_eta_b", positive=True
    )
    k_c = sp.sqrt((alpha - k_tilde + nu_bar * (1 + xi_f0_over_gamma)) / (D_tilde + nu_bar * (1 + xi_f0_over_gamma)))
    k_2 = sp.sqrt((alpha - (1 + eta_s_over_eta_b) * k_tilde) / ((1 + eta_s_over_eta_b) * D_tilde))
    return {
        "k_tilde_trace_zero": k_c,
        "k_tilde_det_root1": sp.Integer(1),
        "k_tilde_det_root2": k_2,
    }
