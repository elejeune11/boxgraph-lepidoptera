import math


def run(inputs, settings):
    a = float(inputs["a"])
    b = float(inputs["b"])
    phi = float(inputs["phi"])
    f = 0.5 * a * phi**2 + 0.25 * b * phi**4
    mu_eq = a * phi + b * phi**3
    P_eq = phi * mu_eq - f
    phi_b = math.sqrt(-a / b)
    phi_s = math.sqrt(-a / (3.0 * b))
    return {
        "f": f,
        "mu_eq": mu_eq,
        "P_eq": P_eq,
        "phi_b_pos": phi_b,
        "phi_b_neg": -phi_b,
        "phi_s_pos": phi_s,
    }


def symbolic():
    import sympy as sp

    a, b, phi = sp.symbols("a b phi", real=True)
    f = sp.Rational(1, 2) * a * phi**2 + sp.Rational(1, 4) * b * phi**4
    mu_eq = a * phi + b * phi**3
    P_eq = sp.simplify(phi * mu_eq - f)
    phi_b = sp.sqrt(-a / b)
    phi_s = sp.sqrt(-a / (3 * b))
    return {
        "f": f,
        "mu_eq": mu_eq,
        "P_eq": P_eq,
        "phi_b_pos": phi_b,
        "phi_b_neg": -phi_b,
        "phi_s_pos": phi_s,
    }
