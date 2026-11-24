def run(inputs, settings):
    # gamma_I in pN, L in um (1000 nm) -> pN*nm
    return {"Delta_E_surface": inputs["gamma_I"] * inputs["L"] * 1000.0}


def symbolic():
    import sympy as sp

    gamma_I, L = sp.symbols("gamma_I L", positive=True)
    return {"Delta_E_surface": 1000 * gamma_I * L}
