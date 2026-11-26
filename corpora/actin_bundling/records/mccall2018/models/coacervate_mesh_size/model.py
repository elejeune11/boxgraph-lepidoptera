def run(inputs, settings):
    return {"xi": inputs["a"] * (inputs["Phi_pp"] ** (-0.75))}


def symbolic():
    import sympy as sp

    a, Phi_pp = sp.symbols("a Phi_pp", positive=True)
    return {"xi": a * Phi_pp ** sp.Rational(-3, 4)}
