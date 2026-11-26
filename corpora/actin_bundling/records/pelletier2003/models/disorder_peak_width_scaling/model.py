def run(inputs, settings):
    # Scaling law printed without a prefactor (prefactor: unknown); unit prefactor = 1 nm.
    q = inputs["q"]
    return {"Delta_q": q * q}


def symbolic():
    import sympy as sp

    q = sp.Symbol("q", positive=True)
    return {"Delta_q": q**2}
