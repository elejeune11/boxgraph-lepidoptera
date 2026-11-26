"""Critical buckling load scaling (p. 7 col 1): F_crit ~ kappa_B / L^2."""


def run(inputs, settings):
    return {"F_crit": inputs["kappa_B"] / (inputs["L"] ** 2)}


def symbolic():
    import sympy as sp

    kappa_B, L = sp.symbols("kappa_B L", positive=True)
    return {"F_crit": kappa_B / L**2}
