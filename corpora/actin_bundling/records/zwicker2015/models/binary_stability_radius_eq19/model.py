def run(inputs, settings):
    return {"R_stab": (3.0 * inputs["D_B"] * inputs["gamma_beta"] / (2.0 * inputs["k_b"])) ** (1.0 / 3.0)}


def symbolic():
    import sympy as sp

    D_B, gamma_beta, k_b = sp.symbols("D_B gamma_beta k_b", positive=True)
    return {"R_stab": (3 * D_B * gamma_beta / (2 * k_b)) ** sp.Rational(1, 3)}
