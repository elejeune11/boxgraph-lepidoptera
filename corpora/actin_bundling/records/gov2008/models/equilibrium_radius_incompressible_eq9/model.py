import math


def run(inputs, settings):
    a = inputs["a"]
    L_p = inputs["L_p"]
    alpha_prime = inputs["alpha_prime"]
    rho = inputs["rho"]
    r_over_a = (math.pi * ((L_p / a) ** 4) / (alpha_prime * rho)) ** (1.0 / 3.0)
    return {"R_0_over_a": r_over_a, "R_0": r_over_a * a}


def symbolic():
    import sympy as sp

    a, L_p, alpha_prime, rho = sp.symbols("a L_p alpha_prime rho", positive=True)
    r_over_a = (sp.pi * (L_p / a) ** 4 / (alpha_prime * rho)) ** sp.Rational(1, 3)
    return {"R_0_over_a": r_over_a, "R_0": r_over_a * a}
