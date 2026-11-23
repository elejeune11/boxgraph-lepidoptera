def run(inputs, settings):
    # prefactor: unknown -> unit prefactor (1) per BoxGraph convention
    D_si = inputs["D"] * 1e-12
    return {"zeta_c": inputs["eta"] / (D_si * inputs["tau"])}


def symbolic():
    import sympy as sp

    eta, D, tau = sp.symbols("eta D tau", positive=True)
    return {"zeta_c": eta / (D * sp.Rational(1, 10**12) * tau)}
