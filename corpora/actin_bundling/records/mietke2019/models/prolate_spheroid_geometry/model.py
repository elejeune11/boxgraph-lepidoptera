import math


def run(inputs, settings):
    e = inputs["e"]
    v_ratio = inputs["V_over_R03"]
    scale = (3.0 * v_ratio / (4.0 * math.pi)) ** (1.0 / 3.0)
    one_minus_e2 = 1.0 - e * e
    a = scale * (one_minus_e2 ** (1.0 / 6.0))
    b = scale * (one_minus_e2 ** (-1.0 / 3.0))
    return {
        "a_over_R0": a,
        "b_over_R0": b,
        "Css_pole_dimless": b / (a * a),
        "Css_eq_dimless": a / (b * b),
        "Cphi_eq_dimless": 1.0 / a,
    }


def symbolic():
    import sympy as sp

    e, V_over_R03 = sp.symbols("e V_over_R03", nonnegative=True)
    scale = (3 * V_over_R03 / (4 * sp.pi)) ** sp.Rational(1, 3)
    a = scale * (1 - e**2) ** sp.Rational(1, 6)
    b = scale * (1 - e**2) ** (-sp.Rational(1, 3))
    return {
        "a_over_R0": a,
        "b_over_R0": b,
        "Css_pole_dimless": b / a**2,
        "Css_eq_dimless": a / b**2,
        "Cphi_eq_dimless": 1 / a,
    }
