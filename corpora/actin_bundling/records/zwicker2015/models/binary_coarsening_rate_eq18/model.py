def run(inputs, settings):
    return {
        "lam": inputs["D_B"] * inputs["gamma_beta"] / (inputs["R_bar"] ** 3)
        - (2.0 / 3.0) * inputs["k_b"]
    }


def symbolic():
    import sympy as sp

    D_B, gamma_beta, R_bar, k_b = sp.symbols("D_B gamma_beta R_bar k_b", positive=True)
    return {"lam": D_B * gamma_beta / R_bar**3 - sp.Rational(2, 3) * k_b}
