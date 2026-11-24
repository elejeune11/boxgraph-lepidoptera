import math


def run(inputs, settings):
    V_bar = inputs["k_f"] * inputs["V_s"] / ((inputs["k_f"] + inputs["k_b"]) * inputs["N"])
    R_bar = (3.0 * V_bar / (4.0 * math.pi)) ** (1.0 / 3.0)
    return {"V_bar": V_bar, "R_bar": R_bar}


def symbolic():
    import sympy as sp

    k_f, k_b, V_s, N = sp.symbols("k_f k_b V_s N", positive=True)
    V_bar = k_f * V_s / ((k_f + k_b) * N)
    R_bar = (3 * V_bar / (4 * sp.pi)) ** sp.Rational(1, 3)
    return {"V_bar": V_bar, "R_bar": R_bar}
