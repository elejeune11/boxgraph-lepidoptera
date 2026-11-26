import math


def run(inputs, settings):
    n_a = float(inputs["n_a"])
    phi = float(inputs["phi"])
    alpha = float(inputs["alpha"])
    eps_kBT = float(inputs["eps_kBT"])
    return {"N_peak": (n_a ** 2) * math.exp(2.0 * alpha * phi * eps_kBT)}


def symbolic():
    import sympy as sp

    n_a, phi, alpha, eps_kBT = sp.symbols("n_a phi alpha eps_kBT", positive=True)
    return {"N_peak": n_a**2 * sp.exp(2 * alpha * phi * eps_kBT)}
