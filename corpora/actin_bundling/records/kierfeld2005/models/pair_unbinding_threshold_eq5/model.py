def run(inputs, settings):
    T = inputs["T"]
    L_p = inputs["L_p"]
    l_a = inputs["l_a"]
    return {"W_bar_c": (T / L_p) * (L_p / l_a) ** (2.0 / 3.0)}


def symbolic():
    import sympy as sp
    T, L_p, l_a = sp.symbols("T L_p l_a", positive=True)
    return {"W_bar_c": (T / L_p) * (L_p / l_a) ** sp.Rational(2, 3)}
