import math


def run(inputs, settings):
    V_ratio = inputs["V_ratio"]
    x_s = (0.5 * V_ratio) ** (1.0 / 3.0)
    W_sym = 2.0 * ((1.0 + x_s ** 2) ** 2.5) / (1.0 + 4.0 * x_s ** 2)

    # Unequal-to-single boundary: dF_tot/dR_1 = 0 with smaller droplet radius R_1 = a (V_1 = V_c)
    if abs(V_ratio - 2.0) < 1e-10:
        W_unequal = 8.0 * math.sqrt(2.0) / 5.0
    elif V_ratio > 2.0:
        x2 = (V_ratio - 1.0) ** (1.0 / 3.0)
        num = 2.0 * (1.0 - 1.0 / x2)
        den = 1.0 / (2.0 * math.sqrt(2.0)) - 1.0 / (x2 * ((1.0 + x2 ** 2) ** 1.5))
        W_unequal = num / den
    else:
        W_unequal = W_sym

    return {"W_sym": W_sym, "W_unequal": W_unequal}


def symbolic():
    import sympy as sp
    V_ratio = sp.symbols("V_ratio", positive=True)
    x_s = (V_ratio / 2) ** sp.Rational(1, 3)
    x2 = (V_ratio - 1) ** sp.Rational(1, 3)
    W_sym = 2 * (1 + x_s ** 2) ** sp.Rational(5, 2) / (1 + 4 * x_s ** 2)
    W_unequal = 2 * (1 - 1 / x2) / (1 / (2 * sp.sqrt(2)) - 1 / (x2 * (1 + x2 ** 2) ** sp.Rational(3, 2)))
    return {"W_sym": W_sym, "W_unequal": W_unequal}
