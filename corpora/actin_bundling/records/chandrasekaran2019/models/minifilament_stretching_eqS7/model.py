def run(inputs, settings):
    dl = inputs["l_ij"] - inputs["l_0"]
    return {"U_MF": 0.5 * inputs["K_str_MF"] * dl * dl}


def symbolic():
    import sympy as sp
    K_str_MF, l_ij, l_0 = sp.symbols("K_str_MF l_ij l_0", positive=True)
    return {"U_MF": sp.Rational(1, 2) * K_str_MF * (l_ij - l_0) ** 2}
