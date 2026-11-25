def run(inputs, settings):
    k_tilde = inputs["n_mode"] / inputs["L_over_2pi_r0"]
    eta = inputs["eta_s_over_eta_b"]
    nu_bar = (eta + 1.0) / (4.0 * eta)
    return {"max_re_lambda_tilde": nu_bar * (1.0 - k_tilde * k_tilde)}


def symbolic():
    import sympy as sp

    L_over_2pi_r0, n_mode, eta_s_over_eta_b = sp.symbols(
        "L_over_2pi_r0 n_mode eta_s_over_eta_b", positive=True
    )
    k_tilde = n_mode / L_over_2pi_r0
    return {
        "max_re_lambda_tilde": ((eta_s_over_eta_b + 1) / (4 * eta_s_over_eta_b)) * (1 - k_tilde**2)
    }
