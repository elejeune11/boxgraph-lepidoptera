def run(inputs, settings):
    R = inputs["R"]
    mu_par_bar = inputs["mu_par_bar"]
    gamma = inputs["gamma"]
    theta_0 = gamma / (2.0 * mu_par_bar * R)
    e_elast = -(gamma**2) / (4.0 * mu_par_bar * (R**2))
    return {"theta_0": theta_0, "E_elast": e_elast}


def symbolic():
    import sympy as sp
    R, mu_par_bar, gamma = sp.symbols("R mu_par_bar gamma", positive=True)
    return {
        "theta_0": gamma / (2 * mu_par_bar * R),
        "E_elast": -(gamma**2) / (4 * mu_par_bar * R**2),
    }
