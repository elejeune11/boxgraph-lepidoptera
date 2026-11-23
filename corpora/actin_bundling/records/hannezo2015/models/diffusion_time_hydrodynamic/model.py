def run(inputs, settings):
    # eta [Pa*s*m] / (xi [Pa*s/m] * D [um^2/s -> m^2/s = 1e-12]) -> s
    D_m2_s = inputs["D"] * 1e-12
    return {"T": inputs["eta"] / (inputs["xi"] * D_m2_s)}


def symbolic():
    import sympy as sp

    eta, xi, D = sp.symbols("eta xi D", positive=True)
    return {"T": eta / (xi * D * sp.Rational(1, 10**12))}
