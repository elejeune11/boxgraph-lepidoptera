def run(inputs, settings):
    # Prefactor is unknown (~), so code uses prefactor = 1 as specified by BoxGraph conventions.
    return {"tau": (inputs["mu"] / inputs["kappa"]) * (inputs["L"] ** 4)}


def symbolic():
    import sympy as sp

    mu, kappa, L = sp.symbols("mu kappa L", positive=True)
    return {"tau": (mu / kappa) * L**4}
