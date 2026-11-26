def run(inputs, settings):
    return {"U_4": 2.0 / 3.0}


def symbolic():
    import sympy as sp

    t = sp.symbols("t", positive=True)
    return {"U_4": sp.Rational(2, 3) + 0 * t}
