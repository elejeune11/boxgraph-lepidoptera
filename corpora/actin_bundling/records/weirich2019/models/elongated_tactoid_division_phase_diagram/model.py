def run(inputs, settings):
    V_ratio = inputs["V_ratio"]
    W_el_sym = 1.5 * ((0.5 * V_ratio) ** 0.4)
    if V_ratio >= 2.0:
        a2 = (V_ratio - 1.0) ** 0.2
        W_el_unequal = 3.0 * (a2 ** 2) / (1.0 + a2 ** 2)
    else:
        W_el_unequal = W_el_sym
    return {"W_el_sym": W_el_sym, "W_el_unequal": W_el_unequal}


def symbolic():
    import sympy as sp
    V_ratio = sp.symbols("V_ratio", positive=True)
    a2 = (V_ratio - 1) ** sp.Rational(1, 5)
    return {
        "W_el_sym": sp.Rational(3, 2) * (V_ratio / 2) ** sp.Rational(2, 5),
        "W_el_unequal": 3 * a2 ** 2 / (1 + a2 ** 2),
    }
