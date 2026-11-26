def run(inputs, settings):
    L = inputs["L"]
    return {"rho": 1.0 / (L * L)}


def symbolic():
    import sympy as sp

    L = sp.symbols("L", positive=True)
    return {"rho": 1 / L**2}
