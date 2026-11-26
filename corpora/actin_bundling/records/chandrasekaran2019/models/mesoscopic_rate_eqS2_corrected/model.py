def run(inputs, settings):
    n = settings.get("n", 2)
    return {"c_mu": inputs["k_mu"] / ((inputs["V_r"] * inputs["N_A"]) ** (n - 1))}


def symbolic():
    import sympy as sp
    k_mu, V_r, N_A = sp.symbols("k_mu V_r N_A", positive=True)
    return {"c_mu": k_mu / (V_r * N_A)}
