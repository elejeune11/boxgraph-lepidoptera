def run(inputs, settings):
    W_el = float(inputs["W_el"])
    V_ratio_sym = 2.0 * ((2.0 * W_el / 3.0) ** 2.5)
    return {"V_ratio_sym": V_ratio_sym}


def symbolic():
    import sympy as sp
    W_el = sp.symbols("W_el", positive=True)
    return {"V_ratio_sym": 2 * (2 * W_el / 3) ** sp.Rational(5, 2)}
