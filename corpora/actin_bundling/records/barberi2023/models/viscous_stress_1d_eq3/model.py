def run(inputs, settings):
    eta = float(inputs["eta"])
    eta_tilde = float(inputs["eta_tilde"])
    d = float(inputs["d"])
    dv_dx = float(inputs["dv_dx"])
    # Printed p. 1 formula: v_tensor = (1/2) * [grad v + (grad v)^T - (div v / d) I]
    # In 1D along x: (1/2) * (2 - 1/d) * dv_dx, so 2 * eta * v_tensor = eta * (2 - 1/d) * dv_dx.
    eff_visc = eta * (2.0 - 1.0 / d) + eta_tilde
    return {
        "sigma_visc": eff_visc * dv_dx,
        "viscous_coeff": eff_visc / eta,
    }


def symbolic():
    import sympy as sp

    eta, eta_tilde, d, dv_dx = sp.symbols("eta eta_tilde d dv_dx", positive=True)
    eff_visc = eta * (2 - 1 / d) + eta_tilde
    return {
        "sigma_visc": eff_visc * dv_dx,
        "viscous_coeff": eff_visc / eta,
    }
