def run(inputs, settings):
    dl = inputs["l_ij"] - inputs["l_0"]
    return {"U_linker": 0.5 * inputs["K_str_linker"] * dl * dl}


def symbolic():
    import sympy as sp
    K_str_linker, l_ij, l_0 = sp.symbols("K_str_linker l_ij l_0", positive=True)
    return {"U_linker": sp.Rational(1, 2) * K_str_linker * (l_ij - l_0) ** 2}
