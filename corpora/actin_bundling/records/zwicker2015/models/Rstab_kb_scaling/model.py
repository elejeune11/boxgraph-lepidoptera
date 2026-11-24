def run(inputs, settings):
    return {"R_stab": inputs["k_b"] ** (-1.0 / 3.0)}


def symbolic():
    import sympy as sp

    k_b = sp.symbols("k_b", positive=True)
    return {"R_stab": k_b ** (-sp.Rational(1, 3))}
