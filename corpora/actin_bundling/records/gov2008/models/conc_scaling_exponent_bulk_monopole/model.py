import math
from scipy.special import lambertw


def run(inputs, settings):
    c = inputs["c"]
    eps = inputs["eps"]
    eps_b_max = inputs["eps_b_max"]
    k_b_prime = inputs["k_b_prime"]
    a = inputs["a"]
    rho = inputs["rho"]
    ce = c * math.exp(eps)
    eps_b = eps_b_max * ce / (1.0 + ce)
    u = (a * a) * rho
    C = 3.0 * math.exp(1.5) * eps_b / (2.0 * k_b_prime * u * u)
    w = float(lambertw(C).real)
    return {"exp_c": ((1.0 / 3.0) * w / (1.0 + w)) / (1.0 + ce)}


def symbolic():
    import sympy as sp

    c, eps, eps_b_max, k_b_prime, a, rho = sp.symbols("c eps eps_b_max k_b_prime a rho", positive=True)
    ce = c * sp.exp(eps)
    eps_b = eps_b_max * ce / (1 + ce)
    C = 3 * sp.exp(sp.Rational(3, 2)) * eps_b / (2 * k_b_prime * (a**2 * rho) ** 2)
    W = sp.LambertW(C)
    return {"exp_c": (sp.Rational(1, 3) * W / (1 + W)) / (1 + ce)}
