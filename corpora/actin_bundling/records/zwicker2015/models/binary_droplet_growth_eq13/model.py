def run(inputs, settings):
    R = inputs["R"]
    D_B = inputs["D_B"]
    phi_0_B = inputs["phi_0_B"]
    gb = inputs["gamma_beta"]
    k_b = inputs["k_b"]
    return {"dR_dt": D_B * (phi_0_B / R - gb / R**2) - k_b * R / 3.0}


def symbolic():
    import sympy as sp

    R, D_B, phi_0_B, gamma_beta, k_b = sp.symbols("R D_B phi_0_B gamma_beta k_b", positive=True)
    return {"dR_dt": D_B * (phi_0_B / R - gamma_beta / R**2) - k_b * R / 3}
