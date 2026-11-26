def run(inputs, settings):
    dl = inputs["l_m"] - inputs["l_0m"]
    return {"U_motor": 0.5 * inputs["K_motor"] * dl * dl}

def symbolic():
    import sympy as sp
    K_motor, l_m, l_0m = sp.symbols("K_motor l_m l_0m", real=True)
    return {"U_motor": sp.Rational(1, 2) * K_motor * (l_m - l_0m) ** 2}
