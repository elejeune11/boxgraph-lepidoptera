def run(inputs, settings):
    return {"R_crit": inputs["gamma_beta"] / inputs["phi_0_B"]}


def symbolic():
    import sympy as sp

    gamma_beta, phi_0_B = sp.symbols("gamma_beta phi_0_B", positive=True)
    return {"R_crit": gamma_beta / phi_0_B}
