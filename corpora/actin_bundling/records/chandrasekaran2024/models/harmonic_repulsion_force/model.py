def run(inputs, settings):
    # k_spring in pN/um; delta in nm -> 1e-3 um
    d = max(0.0, inputs["delta"])
    return {"F": inputs["k_spring"] * (d * 1e-3)}


def symbolic():
    import sympy as sp

    k_spring, delta = sp.symbols("k_spring delta", positive=True)
    return {"F": sp.Rational(1, 1000) * k_spring * delta}
