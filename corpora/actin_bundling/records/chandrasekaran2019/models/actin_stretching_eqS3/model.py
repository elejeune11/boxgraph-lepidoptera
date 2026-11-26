def run(inputs, settings):
    dl = inputs["l_i"] - inputs["l_0"]
    return {"U_str": 0.5 * inputs["K_str"] * dl * dl}


def symbolic():
    import sympy as sp
    K_str, l_i, l_0 = sp.symbols("K_str l_i l_0", positive=True)
    return {"U_str": sp.Rational(1, 2) * K_str * (l_i - l_0) ** 2}
