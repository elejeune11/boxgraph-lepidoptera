import math


def run(inputs, settings):
    c = inputs["c"]
    eps = inputs["eps"]
    eps_b_max = inputs["eps_b_max"]
    rho = inputs["rho"]
    a = inputs["a"]
    lam = inputs["lam"]
    E_c = inputs["E_c"]
    ce = c * math.exp(eps)
    eps_b = eps_b_max * ce / (1.0 + ce)
    ec_star = eps_b / (rho * a * lam)
    if E_c <= ec_star:
        return {"exp_c": None}
    ratio = ec_star / E_c
    return {"exp_c": (ratio / (1.0 - ratio)) / (1.0 + ce)}


def symbolic():
    import sympy as sp

    c, eps, eps_b_max, rho, a, lam, E_c = sp.symbols("c eps eps_b_max rho a lam E_c", positive=True)
    ce = c * sp.exp(eps)
    eps_b = eps_b_max * ce / (1 + ce)
    ratio = eps_b / (rho * a * lam * E_c)
    return {"exp_c": sp.simplify((ratio / (1 - ratio)) / (1 + ce))}
