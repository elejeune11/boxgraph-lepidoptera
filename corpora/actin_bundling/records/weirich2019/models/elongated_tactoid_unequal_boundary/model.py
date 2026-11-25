def run(inputs, settings):
    W_el = float(inputs["W_el"])
    V_ratio_unequal = 1.0 + ((W_el / (3.0 - W_el)) ** 2.5)
    return {"V_ratio_unequal": V_ratio_unequal}


def symbolic():
    import sympy as sp
    W_el = sp.symbols("W_el", positive=True)
    return {"V_ratio_unequal": 1 + (W_el / (3 - W_el)) ** sp.Rational(5, 2)}
