def run(inputs, settings):
    r = inputs["r"]
    r_c = inputs["r_c"]
    if r >= r_c:
        u = 0.5 * inputs["k_c"] * (r - r_c) ** 2
    else:
        u = 0.0
    return {"U_confine": u}


def symbolic():
    import sympy as sp

    k_c, r, r_c = sp.symbols("k_c r r_c", nonnegative=True)
    return {"U_confine": sp.Piecewise((sp.Rational(1, 2) * k_c * (r - r_c) ** 2, r >= r_c), (0, True))}
