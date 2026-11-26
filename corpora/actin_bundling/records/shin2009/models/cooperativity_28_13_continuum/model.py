def run(inputs, settings):
    beta_U = float(inputs["beta_U"])
    xi_t = float(inputs["xi_t"])
    n0 = float(inputs["f_bound"])
    dw = float(inputs["delta_omega"])
    # Integral: 2 * (1/2) * n0 * int_0^inf cos(dw * x) * exp(-x / xi_t) dx = n0 * xi_t / (1 + (dw * xi_t)^2)
    denom = 1.0 + (dw * xi_t) ** 2
    u_norm_cont = n0 * xi_t / denom
    u_cont = (beta_U ** 2) * u_norm_cont
    xi_t_max = 1.0 / dw
    local_xit_exp = (1.0 - (dw * xi_t) ** 2) / denom
    return {
        "u_cont": u_cont,
        "u_norm_cont": u_norm_cont,
        "xi_t_max": xi_t_max,
        "local_xit_exponent": local_xit_exp,
    }


def symbolic():
    import sympy as sp

    beta_U, xi_t, f_bound, delta_omega = sp.symbols(
        "beta_U xi_t f_bound delta_omega", positive=True
    )
    denom = 1 + (delta_omega * xi_t) ** 2
    u_norm_cont = f_bound * xi_t / denom
    return {
        "u_cont": beta_U ** 2 * u_norm_cont,
        "u_norm_cont": u_norm_cont,
        "xi_t_max": 1 / delta_omega,
        "local_xit_exponent": (1 - (delta_omega * xi_t) ** 2) / denom,
    }
