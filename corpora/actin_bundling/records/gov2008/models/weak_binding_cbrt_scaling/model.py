def run(inputs, settings):
    return {"R_0": inputs["eps_b"] ** (1.0 / 3.0)}


def symbolic():
    import sympy as sp

    eps_b = sp.symbols("eps_b", positive=True)
    return {"R_0": eps_b ** sp.Rational(1, 3)}
