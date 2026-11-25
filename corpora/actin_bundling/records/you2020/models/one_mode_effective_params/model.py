def run(inputs, settings):
    alpha_A = float(inputs["chi_A"]) + float(inputs["gamma_A"]) + float(inputs["phi0_A"]) ** 2
    alpha_B = float(inputs["chi_B"]) + float(inputs["gamma_B"]) + float(inputs["phi0_B"]) ** 2
    return {"alpha_A": alpha_A, "alpha_B": alpha_B}


def symbolic():
    import sympy as sp

    chi_A, gamma_A, phi0_A, chi_B, gamma_B, phi0_B = sp.symbols(
        "chi_A gamma_A phi0_A chi_B gamma_B phi0_B", real=True
    )
    return {
        "alpha_A": chi_A + gamma_A + phi0_A**2,
        "alpha_B": chi_B + gamma_B + phi0_B**2,
    }
