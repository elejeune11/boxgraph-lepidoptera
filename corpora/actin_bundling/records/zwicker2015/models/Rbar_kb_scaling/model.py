def run(inputs, settings):
    return {"R_bar": inputs["k_b"] ** (-0.5)}


def symbolic():
    import sympy as sp

    k_b = sp.symbols("k_b", positive=True)
    return {"R_bar": k_b ** (-sp.Rational(1, 2))}
