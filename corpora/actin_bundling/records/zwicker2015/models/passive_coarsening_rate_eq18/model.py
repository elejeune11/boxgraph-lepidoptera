def run(inputs, settings):
    return {"lam": inputs["D_B"] * inputs["gamma_beta"] / (inputs["R_bar"] ** 3)}


def symbolic():
    import sympy as sp

    D_B, gamma_beta, R_bar = sp.symbols("D_B gamma_beta R_bar", positive=True)
    return {"lam": D_B * gamma_beta / R_bar**3}
