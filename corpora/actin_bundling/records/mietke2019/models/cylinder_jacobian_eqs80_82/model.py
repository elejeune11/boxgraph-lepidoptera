import math


def jacobian_entries(k_tilde, alpha, xi_f0, kappa_tilde, D_tilde, k_turn, eta_ratio):
    k2 = k_tilde * k_tilde
    # Eq. (S79): r0^4 * B_c(k_tilde) = 4*k^2 + (2*k^2 - 3)*(k^2 - 1) = 2*k^4 - k^2 + 3
    bc = 4.0 * k2 + (2.0 * k2 - 3.0) * (k2 - 1.0)
    term = -kappa_tilde * bc + (1.0 + xi_f0) * (1.0 - k2)
    nu_bar = (eta_ratio + 1.0) / (4.0 * eta_ratio)
    J11 = nu_bar * term
    J12 = -0.5 * alpha
    J21 = -0.5 * term
    J22 = alpha - D_tilde * k2 - k_turn
    return J11, J12, J21, J22


def run(inputs, settings):
    k_tilde = inputs["n_mode"] / inputs["L_over_2pi_r0"]
    J11, J12, J21, J22 = jacobian_entries(
        k_tilde,
        inputs["alpha"],
        inputs["xi_f0_over_gamma"],
        inputs["kappa_tilde"],
        inputs["D_tilde"],
        inputs["k_tilde"],
        inputs["eta_s_over_eta_b"],
    )
    tr = J11 + J22
    det = J11 * J22 - J12 * J21
    disc = tr * tr - 4.0 * det
    if disc >= 0.0:
        max_re = 0.5 * (tr + math.sqrt(disc))
        im = 0.0
        period = 0.0
    else:
        max_re = 0.5 * tr
        im = 0.5 * math.sqrt(-disc)
        period = 2.0 * math.pi / im if im > 0.0 else 0.0
    return {
        "max_re_lambda_tilde": max_re,
        "im_lambda_tilde": im,
        "oscillation_period_tilde": period,
        "tr_tilde": tr,
        "det_tilde": det,
    }
