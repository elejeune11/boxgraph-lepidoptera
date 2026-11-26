def run(inputs, settings):
    n = settings.get("n", 2)
    # Printed Eq. (S2): k_mu * (V_r / N_A)^(n - 1).
    # At n=2: [L/(mol*s)] * [L / (1/mol)] = L^2 / s (dimensionally NOT 1/s!).
    return {"c_mu_printed": inputs["k_mu"] * (inputs["V_r"] / inputs["N_A"]) ** (n - 1)}


def symbolic():
    import sympy as sp
    k_mu, V_r, N_A = sp.symbols("k_mu V_r N_A", positive=True)
    return {"c_mu_printed": k_mu * (V_r / N_A)}
