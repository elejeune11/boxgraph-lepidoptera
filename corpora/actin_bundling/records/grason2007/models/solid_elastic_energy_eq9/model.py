def run(inputs, settings):
    R = inputs["R"]
    Omega = inputs["Omega"]
    mu_par_bar = inputs["mu_par_bar"]
    mu_perp_bar = inputs["mu_perp_bar"]
    K3_bar = inputs["K3_bar"]
    gamma = inputs["gamma"]
    E = (
        mu_par_bar * (Omega * R) ** 2
        + mu_perp_bar * (Omega * R) ** 4
        + K3_bar * (Omega**4) * (R**2)
        - gamma * Omega
    )
    return {"E_elast": E}


def symbolic():
    import sympy as sp
    R, Omega, mu_par_bar, mu_perp_bar, K3_bar, gamma = sp.symbols(
        "R Omega mu_par_bar mu_perp_bar K3_bar gamma", positive=True
    )
    return {
        "E_elast": mu_par_bar * (Omega * R) ** 2
        + mu_perp_bar * (Omega * R) ** 4
        + K3_bar * Omega**4 * R**2
        - gamma * Omega
    }
