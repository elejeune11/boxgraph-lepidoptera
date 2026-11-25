def run(inputs, settings):
    # D [um^2/s] * beta [Pa*s] / L [pN = Pa*um^2] = dimensionless
    return {"delta": inputs["D"] * inputs["kappa"] * inputs["beta"] / (2.0 * inputs["L"])}

def symbolic():
    import sympy as sp
    D, kappa, beta, L = sp.symbols("D kappa beta L", real=True)
    return {"delta": D * kappa * beta / (2 * L)}
