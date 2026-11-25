import math


def _solve_mode(inputs, use_corrected_diffusion):
    l = int(round(inputs["l"]))
    eta_s = inputs["eta_s"]
    eta_b = inputs["eta_b"]
    gamma = inputs["gamma"]
    xi = inputs["xi"]
    f0 = inputs["f_c0"]
    c0_df = inputs["c0_df_dc0"]
    kappa = inputs["kappa"]
    R0 = inputs["R0"]
    D = inputs["D"]
    k = inputs["k"]

    # Eq. (S71) as printed has [R0^2 * D * l*(l+1) + k]; corrected has [D * l*(l+1) / R0^2 + k].
    if use_corrected_diffusion:
        diff_term = D * l * (l + 1) / (R0 ** 2)
    else:
        diff_term = (R0 ** 2) * D * l * (l + 1)

    if l == 1:
        lam1 = -diff_term - k + (xi / eta_b) * c0_df
        return {"max_re_lambda": lam1, "im_lambda": 0.0}

    # For l >= 2, exact 2x2 Jacobian in (delta R_l / R0, delta c_l / c0) from Eqs. (S69)-(S71):
    # Effective mode stiffness S_l = gamma + xi*f(c0) + (kappa/R0^2)*l*(l+1)
    S_l = gamma + xi * f0 + (kappa / (R0 ** 2)) * l * (l + 1)
    J11 = -S_l * (eta_b * l * (l + 1) + eta_s * (l - 1) * (l + 2)) / (4.0 * eta_b * eta_s)
    J12 = -xi * c0_df / (2.0 * eta_b)
    J21 = (l - 1) * (l + 2) * S_l / (2.0 * eta_b)
    J22 = (xi * c0_df / eta_b) - diff_term - k

    tr = J11 + J22
    det = J11 * J22 - J12 * J21
    disc = tr * tr - 4.0 * det
    if disc >= 0.0:
        max_re = 0.5 * (tr + math.sqrt(disc))
        im = 0.0
    else:
        max_re = 0.5 * tr
        im = 0.5 * math.sqrt(-disc)
    return {"max_re_lambda": max_re, "im_lambda": im}


def run(inputs, settings):
    return _solve_mode(inputs, use_corrected_diffusion=False)
