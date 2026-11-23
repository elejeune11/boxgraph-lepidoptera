def run(inputs, settings):
    R = inputs["R"]
    Omega = inputs["Omega"]
    mu_perp_bar = inputs["mu_perp_bar"]
    K2_bar = inputs["K2_bar"]
    K3_bar = inputs["K3_bar"]
    gamma = inputs["gamma"]
    E = (
        mu_perp_bar * (Omega * R) ** 4
        + K2_bar * Omega**2
        + K3_bar * (Omega**4) * (R**2)
        - gamma * Omega
    )
    return {"E_elast": E}


def symbolic():
    import sympy as sp
    R, Omega, mu_perp_bar, K2_bar, K3_bar, gamma = sp.symbols(
        "R Omega mu_perp_bar K2_bar K3_bar gamma", positive=True
    )
    return {
        "E_elast": mu_perp_bar * (Omega * R) ** 4
        + K2_bar * Omega**2
        + K3_bar * Omega**4 * R**2
        - gamma * Omega
    }
