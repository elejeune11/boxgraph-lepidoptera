def run(inputs, settings):
    W_el = inputs["W_el"]
    V_ratio = inputs["V_ratio"]
    f_V = inputs["f_V"]
    K = inputs["K"]                       # pN
    gamma_pN_um = inputs["gamma"] * 1e-3  # nN/m -> pN/um

    a1 = (max(0.0, f_V) * V_ratio) ** 0.2
    a2 = (max(0.0, 1.0 - f_V) * V_ratio) ** 0.2
    F_tot_dimless = (a1 ** 3) + (a2 ** 3) - W_el * (a1 + a2)
    F_tot = ((K ** 2) / gamma_pN_um) * F_tot_dimless

    # Symmetric second derivative d^2(F_tot/(K^2/gamma))/d(alpha_1)^2 at alpha_1 = alpha_2 = alpha*
    a_s = (0.5 * V_ratio) ** 0.2
    d2F_sym = -12.0 * a_s + 8.0 * W_el / a_s
    return {
        "F_tot_dimless": F_tot_dimless,
        "F_tot": F_tot,
        "d2F_sym": d2F_sym,
    }


def symbolic():
    import sympy as sp
    W_el, V_ratio, f_V, K, gamma = sp.symbols("W_el V_ratio f_V K gamma", positive=True)
    a1 = (f_V * V_ratio) ** sp.Rational(1, 5)
    a2 = ((1 - f_V) * V_ratio) ** sp.Rational(1, 5)
    F_tot_dimless = a1 ** 3 + a2 ** 3 - W_el * (a1 + a2)
    a_s = (V_ratio / 2) ** sp.Rational(1, 5)
    d2F_sym = -12 * a_s + 8 * W_el / a_s
    return {
        "F_tot_dimless": F_tot_dimless,
        "F_tot": (K ** 2 / (gamma * sp.Rational(1, 1000))) * F_tot_dimless,
        "d2F_sym": d2F_sym,
    }
