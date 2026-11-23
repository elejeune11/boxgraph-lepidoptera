def run(inputs, settings):
    # eta [Pa*s*m] / (D [um^2/s -> 1e-12 m^2/s] * tau [s] * xi [Pa*s/m]) -> dimensionless
    return {"phi": inputs["eta"] / (inputs["D"] * 1e-12 * inputs["tau"] * inputs["xi"])}


def symbolic():
    import sympy as sp

    eta, D, tau, xi = sp.symbols("eta D tau xi", positive=True)
    return {"phi": eta / (D * sp.Rational(1, 10**12) * tau * xi)}
