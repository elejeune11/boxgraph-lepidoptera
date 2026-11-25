def run(inputs, settings):
    # gamma [Pa*s/m^2] * D_si [m^2/s] * beta [Pa*s] / (eta [Pa*s] * L_si [Pa*m^2]) -> [1/m^2]
    # Multiply by 1e-12 to express in [1/um^2].
    D_si = inputs["D"] * 1e-12
    L_si = inputs["L"] * 1e-12
    val_m2 = inputs["gamma"] * D_si * inputs["kappa"] * inputs["beta"] / (2.0 * inputs["eta"] * L_si)
    return {"delta_printed": val_m2 * 1e-12}

def symbolic():
    import sympy as sp
    gamma, D, kappa, beta, eta, L = sp.symbols("gamma D kappa beta eta L", real=True)
    return {"delta_printed": gamma * D * kappa * beta / (2 * eta * L) * sp.Rational(1, 10**12)}
