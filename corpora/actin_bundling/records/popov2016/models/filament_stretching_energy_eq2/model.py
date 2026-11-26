def run(inputs, settings):
    dl = inputs["l"] - inputs["l_0"]
    return {"U_str": 0.5 * inputs["K_str"] * dl * dl}

def symbolic():
    import sympy as sp
    K_str, l, l_0 = sp.symbols("K_str l l_0", real=True)
    return {"U_str": sp.Rational(1, 2) * K_str * (l - l_0) ** 2}
