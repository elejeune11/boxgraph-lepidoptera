def run(inputs, settings):
    return {"R_star": inputs["n_f"] ** (5.0 / 3.0)}


def symbolic():
    import sympy as sp

    n_f = sp.Symbol("n_f", positive=True)
    return {"R_star": n_f ** sp.Rational(5, 3)}
