def run(inputs, settings):
    phi_0_B = inputs["phibar"] - inputs["k_b"] / inputs["k"]
    if phi_0_B <= 0.0:
        return {"R_crit": None}
    return {"R_crit": inputs["gamma_beta"] / phi_0_B}


def symbolic():
    import sympy as sp

    gamma_beta, phibar, k_b, k = sp.symbols("gamma_beta phibar k_b k", positive=True)
    return {"R_crit": gamma_beta / (phibar - k_b / k)}
