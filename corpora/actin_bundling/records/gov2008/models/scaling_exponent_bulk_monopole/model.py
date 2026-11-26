import math
from scipy.special import lambertw


def run(inputs, settings):
    eps_b = inputs["eps_b"]
    k_b_prime = inputs["k_b_prime"]
    a = inputs["a"]
    rho = inputs["rho"]
    u = (a * a) * rho
    alpha_b = k_b_prime / eps_b
    C = 3.0 * math.exp(1.5) / (2.0 * alpha_b * u * u)
    w = float(lambertw(C).real)
    return {"exp_eps_b": (1.0 / 3.0) * w / (1.0 + w)}


def symbolic():
    import sympy as sp

    eps_b, k_b_prime, a, rho = sp.symbols("eps_b k_b_prime a rho", positive=True)
    C = 3 * sp.exp(sp.Rational(3, 2)) * eps_b / (2 * k_b_prime * (a**2 * rho) ** 2)
    W = sp.LambertW(C)
    return {"exp_eps_b": sp.Rational(1, 3) * W / (1 + W)}
