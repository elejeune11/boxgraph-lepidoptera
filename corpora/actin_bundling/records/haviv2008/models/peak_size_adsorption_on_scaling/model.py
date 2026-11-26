import math


def run(inputs, settings):
    n_a = float(inputs["n_a"])
    phi = float(inputs["phi"])
    beta = float(inputs["beta"])
    alpha = float(inputs["alpha"])
    eps_kBT = float(inputs["eps_kBT"])
    w = beta * phi * math.exp(eps_kBT)
    p_fa = w / (1.0 + w)
    return {"N_peak": ((n_a * p_fa) ** 2) * math.exp(2.0 * alpha * phi * eps_kBT)}


def symbolic():
    import sympy as sp

    n_a, phi, beta, alpha, eps_kBT = sp.symbols(
        "n_a phi beta alpha eps_kBT", positive=True
    )
    w = beta * phi * sp.exp(eps_kBT)
    return {"N_peak": (n_a * w / (1 + w)) ** 2 * sp.exp(2 * alpha * phi * eps_kBT)}
