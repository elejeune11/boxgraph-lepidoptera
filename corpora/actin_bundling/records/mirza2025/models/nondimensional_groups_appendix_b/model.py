def run(inputs, settings):
    eta = inputs["eta"]
    ell_s = inputs["ell_s"]
    D = inputs["D"]
    k_d = inputs["k_d"]
    er = inputs["eta_rot_ratio"]
    br = inputs["beta_ratio"]
    ar = inputs["a_over_eta_rot"]
    b_rate = inputs["b_over_eta_rot"]
    lo_bar = inputs["lambda_o_bar"]
    L_pN = inputs["L"]

    gamma = eta / ((ell_s * 1e-6) ** 2)
    t_0 = (ell_s ** 2) / D
    k_bar_d = k_d * t_0
    a_bar = ar * er * t_0
    b_bar = b_rate * er * t_0
    c_bar_0 = 2.0 * a_bar - lo_bar
    # L [pN = Pa*um^2] / (eta [Pa*s] * D [um^2/s]) = dimensionless
    L_bar = L_pN / (eta * D)
    eta_rot = er * eta
    beta = br * eta
    a = ar * eta_rot
    b = b_rate * eta_rot
    gamma_D = eta / t_0  # Pa
    rho_lam_o = lo_bar * gamma_D
    c_0 = 2.0 * a - rho_lam_o
    return {
        "gamma": gamma,
        "t_0": t_0,
        "k_bar_d": k_bar_d,
        "a_bar": a_bar,
        "b_bar": b_bar,
        "c_bar_0": c_bar_0,
        "L_bar": L_bar,
        "eta_rot": eta_rot,
        "beta": beta,
        "a": a,
        "b": b,
        "c_0": c_0,
        "rho_lam_o": rho_lam_o,
    }

def symbolic():
    import sympy as sp
    eta, ell_s, D, k_d, eta_rot_ratio, beta_ratio, a_over_eta_rot, b_over_eta_rot, lambda_o_bar, L = sp.symbols(
        "eta ell_s D k_d eta_rot_ratio beta_ratio a_over_eta_rot b_over_eta_rot lambda_o_bar L", real=True
    )
    gamma = eta / (ell_s * sp.Rational(1, 10**6)) ** 2
    t_0 = ell_s**2 / D
    k_bar_d = k_d * t_0
    a_bar = a_over_eta_rot * eta_rot_ratio * t_0
    b_bar = b_over_eta_rot * eta_rot_ratio * t_0
    c_bar_0 = 2 * a_bar - lambda_o_bar
    L_bar = L / (eta * D)
    eta_rot = eta_rot_ratio * eta
    beta = beta_ratio * eta
    a = a_over_eta_rot * eta_rot
    b = b_over_eta_rot * eta_rot
    rho_lam_o = lambda_o_bar * (eta / t_0)
    c_0 = 2 * a - rho_lam_o
    return {
        "gamma": gamma,
        "t_0": t_0,
        "k_bar_d": k_bar_d,
        "a_bar": a_bar,
        "b_bar": b_bar,
        "c_bar_0": c_bar_0,
        "L_bar": L_bar,
        "eta_rot": eta_rot,
        "beta": beta,
        "a": a,
        "b": b,
        "c_0": c_0,
        "rho_lam_o": rho_lam_o,
    }
