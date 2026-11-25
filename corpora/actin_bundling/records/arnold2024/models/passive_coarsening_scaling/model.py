def run(inputs, settings):
    # Prefactor is unknown (scaling law); prefactor set to 1 per BoxGraph rules
    return {"a": inputs["t"] ** (1.0 / 3.0)}


def symbolic():
    import sympy as sp

    t = sp.symbols("t", positive=True)
    return {"a": t ** sp.Rational(1, 3)}
