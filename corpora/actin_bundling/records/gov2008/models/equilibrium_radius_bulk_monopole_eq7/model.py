import math
from scipy.special import lambertw


def run(inputs, settings):
    alpha_b = inputs["alpha_b"]
    a = inputs["a"]
    rho = inputs["rho"]
    u = (a * a) * rho
    denom_coeff = 2.0 * alpha_b * (u * u)
    arg = 3.0 * math.exp(1.5) / denom_coeff
    w = float(lambertw(arg).real)
    r_over_a = (3.0 / (denom_coeff * w)) ** (1.0 / 3.0)
    return {"R_0_over_a": r_over_a, "R_0": r_over_a * a}


def symbolic():
    import sympy as sp

    alpha_b, a, rho = sp.symbols("alpha_b a rho", positive=True)
    u = a**2 * rho
    denom_coeff = 2 * alpha_b * u**2
    arg = 3 * sp.exp(sp.Rational(3, 2)) / denom_coeff
    r_over_a = (3 / (denom_coeff * sp.LambertW(arg))) ** sp.Rational(1, 3)
    return {"R_0_over_a": r_over_a, "R_0": r_over_a * a}
