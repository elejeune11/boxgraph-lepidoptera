def run(inputs, settings):
    shift = inputs["gamma_beta"] / inputs["R"]
    return {
        "phi_plus_B": shift,
        "phi_minus_B": inputs["psi"] + shift,
    }


def symbolic():
    import sympy as sp

    psi, gamma_beta, R = sp.symbols("psi gamma_beta R", positive=True)
    return {
        "phi_plus_B": gamma_beta / R,
        "phi_minus_B": psi + gamma_beta / R,
    }
